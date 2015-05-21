/**
 * $Id: generate.c 1246 2014-06-02 09:07am pdorazio $
 *
 * @brief Red Pitaya simple signal/function generator with pre-defined
 *        signal types.
 *
 * @Author Ales Bardorfer <ales.bardorfer@redpitaya.com>
 *
 * (c) Red Pitaya  http://www.redpitaya.com
 *
 * This part of code is written in C programming language.
 * Please visit http://en.wikipedia.org/wiki/C_(programming_language)
 * for more details on the language used herein.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fpga_awg.h"
#include "version.h"
#include "kiss_fftr.h"

/**
 * GENERAL DESCRIPTION:
 *
 * The code below performs a function of a signal generator, which produces
 * a a signal of user-selectable pred-defined Signal shape
 * [Sine, Square, Triangle], Amplitude and Frequency on a selected Channel:
 *
 *
 *                   /-----\
 *   Signal shape -->|     | -->[data]--+-->[FPGA buf 1]--><DAC 1>
 *   Amplitude ----->| AWG |            |
 *   Frequency ----->|     |             -->[FPGA buf 2]--><DAC 2>
 *                   \-----/            ^
 *                                      |
 *   Channel ---------------------------+ 
 *
 *
 * This is achieved by first parsing the four parameters defining the 
 * signal properties from the command line, followed by synthesizing the 
 * signal in data[] buffer @ 125 MHz sample rate within the
 * generate_signal() function, depending on the Signal shape, Amplitude
 * and Frequency parameters. The data[] buffer is then transferred
 * to the specific FPGA buffer, defined by the Channel parameter -
 * within the write_signal_fpga() function.
 * The FPGA logic repeatably sends the data from both FPGA buffers to the
 * corresponding DACs @ 125 MHz, which in turn produces the synthesized
 * signal on Red Pitaya SMA output connectors labeled DAC1 & DAC2.
 *
 */

/** Maximal signal frequency [Hz] */
const double c_max_frequency = 62.5e6;

/** Minimal signal frequency [Hz] */
const double c_min_frequency = 0;

/** Maximal signal amplitude [Vpp] */
const double c_max_amplitude = 2.0;

/** AWG buffer length [samples]*/
#define n (16*1024)

/** AWG data buffer */
int32_t data[n];

/** Program name */
const char *g_argv0 = NULL;

/** Signal types */
typedef enum {
    eSignalSine,         ///< Sinusoidal waveform.
    eSignalSquare,       ///< Square waveform.
    eSignalTriangle,     ///< Triangular waveform.
    eSignalSweep         ///< Sinusoidal frequency sweep.
} signal_e;

/** AWG FPGA parameters */
typedef struct {
    int32_t  offsgain;   ///< AWG offset & gain.
    uint32_t wrap;       ///< AWG buffer wrap value.
    uint32_t step;       ///< AWG step interval.
} awg_param_t;

/* Synthesizer options */

struct synth_sin_params {
};

struct synth_sqr_params {
    double dcycle;	///< Duty cycle
};

struct synth_tri_params {
};

struct generate_params {
    double ampl;	///< P2P amplitude
    double freq;	///< Frequency
    double dc_off;	///< DC Offset
    double endfreq;	///< End frequency for sweep
    signal_e type;	///< Signal type
    double bw;		///< Signal bandwidth
    union {
        struct synth_sin_params sin;
        struct synth_sqr_params sqr;
        struct synth_tri_params tri;
    } opts;
};

/* Forward declarations */
void synthesize_signal(struct generate_params *gparams, int32_t *data,
                     awg_param_t *awg);
void write_data_fpga(uint32_t ch,
                     const int32_t *data,
                     const awg_param_t *awg);

/** Print usage information */
void usage() {

    const char *format =
        "%s version %s-%s\n"
        "\n"
        "Usage: %s   channel amplitude frequency <type> <end frequency>\n"
        "\n"
        "\tchannel     Channel to generate signal on [1, 2].\n"
        "\tamplitude   Peak-to-peak signal amplitude in Vpp [0.0 - %1.1f].\n"
        "\tfrequency   Signal frequency in Hz [%2.1f - %2.1e].\n"
        "\ttype        Signal type [sine, sqr, tri, sweep].\n"
        "\tend frequency   Sweep-to frequency in Hz [%2.1f - %2.1e].\n"
        "\n";

    fprintf( stderr, format, g_argv0, VERSION_STR, REVISION_STR,
             g_argv0, c_max_amplitude, c_min_frequency, c_max_frequency);
}

static int
parse_nvpair(char *nvpair, char **name, char **value)
{
    char *p;

    p = strchr(nvpair, '=');
    if (p == NULL) {
        return (-1);
    }
    *p = '\0';
    *name = nvpair;
    *value = p + 1;
    return (0);
}

static int
parse_nvparam(struct generate_params *gpp, const char *name, const char *value)
{

    /* Signal type specific params */
    switch (gpp->type) {
    case eSignalSquare:
        if (strcmp(name, "dcycle") == 0) {
            double dcycle;

            dcycle = strtod(value, NULL);
            if (dcycle <= 0.0 || dcycle >= 1.0)
                return (-1);
            gpp->opts.sqr.dcycle = dcycle;
            return (0);
        }
        break;

    default:
        break;
    }
    /* Global params */
    if (strcmp(name, "dcoff") == 0) {
        double dc_off;

        dc_off = strtod(value, NULL);
        if (fabs(dc_off) + gpp->ampl > c_max_amplitude) {
            return (-1);
        }
        gpp->dc_off = dc_off;
        return (0);
    }
    if (strcmp(name, "bw") == 0) {
        double bw;

        bw = strtod(value, NULL);
        if (bw <= 0 || bw > c_max_frequency) {
            return (-1);
        }
        gpp->bw = bw;
        return (0);
    }

    return (-1);
}

/** Signal generator main */
int main(int argc, char *argv[])
{
    int i;
    char *name, *value;
    struct generate_params gparams;

    g_argv0 = argv[0];    

    if ( argc < 4 ) {

        usage();
        return -1;
    }

    /* Channel argument parsing */
    uint32_t ch = atoi(argv[1]) - 1; /* Zero based internally */
    if (ch > 1) {
        fprintf(stderr, "Invalid channel: %s\n", argv[1]);
        usage();
        return -1;
    }

    /* Signal amplitude argument parsing */
    gparams.ampl = strtod(argv[2], NULL);
    if ( (gparams.ampl < 0.0) || (gparams.ampl > c_max_amplitude) ) {
        fprintf(stderr, "Invalid amplitude: %s\n", argv[2]);
        usage();
        return -1;
    }

    /* Signal frequency argument parsing */
    gparams.freq = strtod(argv[3], NULL);
    gparams.endfreq = 0;

    /* Signal type argument parsing */
    gparams.type = eSignalSine;
    gparams.dc_off = 0.0;
    gparams.bw = c_max_frequency;
    if (argc > 4) {
        if ( strcmp(argv[4], "sine") == 0) {
            gparams.type = eSignalSine;
        } else if ( strcmp(argv[4], "sqr") == 0) {
            gparams.type = eSignalSquare;
            gparams.opts.sqr.dcycle = 0.5;
        } else if ( strcmp(argv[4], "tri") == 0) {
            gparams.type = eSignalTriangle;
        } else if ( strcmp(argv[4], "sweep") == 0) {
            gparams.type = eSignalSweep;   
        } else {
            fprintf(stderr, "Invalid signal type: %s\n", argv[4]);
            usage();
            return -1;
        }
    }

    if (argc > 5) {
        for (i = 5; i < argc; i++) {
            if (parse_nvpair(argv[i], &name, &value) == -1)
                break;
            if (parse_nvparam(&gparams, name, value) == -1) {
                fprintf(stderr, "Invalid name/value parameter: %s=%s\n", 
                    name, value);
                usage();
                return (-1);
            }
        }
        if (i < argc) {
            gparams.endfreq = strtod(argv[i], NULL);
        }
    }

    /* Check frequency limits */
    if ( (gparams.freq < c_min_frequency) || (gparams.freq > c_max_frequency ) ) {
        fprintf(stderr, "Invalid frequency: %s\n", argv[3]);
        usage();
        return -1;
    }

    /* Prepare data buffer (calculate from input arguments) */
    awg_param_t awg_params;

    synthesize_signal(&gparams, data, &awg_params);

    /* Write the data to the FPGA and set FPGA AWG state machine */
    write_data_fpga(ch, data, &awg_params);
}

static void
low_pass_filter(double *ddata, double bwlim, double srate)
{
    kiss_fftr_cfg config;
    kiss_fft_cpx *spectrum;
    int i, splen;
    double freq, freq1;

    if (srate / 2.0 <= bwlim) {
        /*
         * Nyquist frequency is already lower than our target
         * limit, nothing to do really.
         */
        return;
    }

    config = kiss_fftr_alloc(n, 0, NULL, NULL);
    splen = (n / 2) + 1;
    spectrum = (kiss_fft_cpx *)malloc(sizeof(kiss_fft_cpx) * splen);
    memset(spectrum, '\0', sizeof(kiss_fft_cpx) * splen);

    kiss_fftr(config, (kiss_fft_scalar *)ddata, spectrum);
    kiss_fftr_free(config);
    freq1 = srate / (double)n;

    for (i = 0; i < splen; i++) {
        /*
         * Scale down spectrum, seems to be necessary for the
         * kiss_fftri to produce correct result.
         */ 
        spectrum[i].r /= (double)n;
        spectrum[i].i /= (double)n;
        freq = freq1 * i;
        if (freq <= bwlim) {
            continue;
        }
        /* Erase any bands above our target cut-off frequency */
        spectrum[i].r = 0.0;
        spectrum[i].i = 0.0;
    }

    config = kiss_fftr_alloc(n, 1, NULL, NULL);
    kiss_fftri(config, spectrum, ddata);
    free(spectrum);
    kiss_fftr_free(config);
    return;
}

/**
 * Synthesize a desired signal.
 *
 * Generates/synthesized  a signal, based on three pre-defined signal
 * types/shapes, signal amplitude & frequency. The data[] vector of 
 * samples at 125 MHz is generated to be re-played by the FPGA AWG module.
 *
 * @param gpp->ampl  Signal amplitude [Vpp].
 * @param gpp->freq  Signal frequency [Hz].
 * @param gpp->type  Signal type/shape [Sine, Square, Triangle].
 * @param data  Returned synthesized AWG data vector.
 * @param awg   Returned AWG parameters.
 *
 */
void synthesize_signal(struct generate_params *gpp, int32_t *data,
                       awg_param_t *awg) {

    uint32_t i;
    double y_thrs, dsample;
    double ddata[n];

    /* Various locally used constants - HW specific parameters */
    const int dcoffs = -155;

    /* This is where frequency is used... */
    awg->offsgain = (dcoffs << 16) + 0x1fff;
    awg->step = round(65536 * gpp->freq/c_awg_smpl_freq * n);
    awg->wrap = round(65536 * (n-1));

    if (gpp->type == eSignalSquare) {
        y_thrs = cos(gpp->opts.sqr.dcycle * M_PI);
    }

    /* Fill data[] with appropriate buffer samples */
    for(i = 0; i < n; i++) {
        
        /* Sine */
        if (gpp->type == eSignalSine) {
            dsample =  0.5 * gpp->ampl * cos(2*M_PI*(double)i/(double)n);
        }
 
        /* Square */
        if (gpp->type == eSignalSquare) {
            dsample = cos(2*M_PI*(double)i/(double)n);
            if (dsample > y_thrs) {
                dsample = 0.5 * gpp->ampl;
            } else {
                dsample = -0.5 * gpp->ampl;
            }
        }
        
        /* Triangle */
        if (gpp->type == eSignalTriangle) {
            dsample = -0.5 * gpp->ampl *(acos(cos(2*M_PI*(double)i/(double)n))/M_PI*2-1);
        }

        /* Sweep */
        /* Loops from i = 0 to n = 16*1024. Generates a sine wave signal that
           changes in frequency as the buffer is filled. */
        double start = 2 * M_PI * gpp->freq;
        double end = 2 * M_PI * gpp->endfreq;
        if (gpp->type == eSignalSweep) {
            double sampFreq = c_awg_smpl_freq; // 125 MHz
            double t = i / sampFreq; // This particular sample
            double T = n / sampFreq; // Wave period = # samples / sample frequency
            /* Actual formula. Frequency changes from start to end. */
            dsample = 0.5 * gpp->ampl * (sin((start*T)/log(end/start) * ((exp(t*log(end/start)/T)-1))));
        }
        /* Add DC offset */
        dsample += gpp->dc_off;
        ddata[i] = dsample;
    }

    /* Apply anti-aliasing/bw filter to the signal */
    low_pass_filter(ddata, gpp->bw, gpp->freq * (double)n);

    /* Convert into format that the DAC can eat */
    for (i = 0; i < n; i++) {
        /* 1 Vpp ==> 8000 DAC counts, from -4000 to 4000 */
        data[i] = round(ddata[i] * 8000.0);
        if (abs(data[i]) > 8191) {
            /* Truncate to max value if needed */
            fprintf(stderr, "Warning: data overflow detected at sample #%d, value %d, truncating to %s8191\n",
              i, data[i], data[i] < 0 ? "-" : "");
            data[i] = data[i] > 0 ? 8191 : -8191;
        }
 
        /* TODO: Remove, not necessary in C/C++. */
        if(data[i] < 0)
            data[i] += (1 << 14);
    }
}

/**
 * Write synthesized data[] to FPGA buffer.
 *
 * @param ch    Channel number [0, 1].
 * @param data  AWG data to write to FPGA.
 * @param awg   AWG paramters to write to FPGA.
 */
void write_data_fpga(uint32_t ch,
                     const int32_t *data,
                     const awg_param_t *awg) {

    uint32_t i;

    fpga_awg_init();

    if(ch == 0) {
        /* Channel A */
        g_awg_reg->state_machine_conf = 0x000041;
        g_awg_reg->cha_scale_off      = awg->offsgain;
        g_awg_reg->cha_count_wrap     = awg->wrap;
        g_awg_reg->cha_count_step     = awg->step;
        g_awg_reg->cha_start_off      = 0;

        for(i = 0; i < n; i++) {
            g_awg_cha_mem[i] = data[i];
        }
    } else {
        /* Channel B */
        g_awg_reg->state_machine_conf = 0x410000;
        g_awg_reg->chb_scale_off      = awg->offsgain;
        g_awg_reg->chb_count_wrap     = awg->wrap;
        g_awg_reg->chb_count_step     = awg->step;
        g_awg_reg->chb_start_off      = 0;

        for(i = 0; i < n; i++) {
            g_awg_chb_mem[i] = data[i];
        }
    }

    /* Enable both channels */
    /* TODO: Should this only happen for the specified channel?
     *       Otherwise, the not-to-be-affected channel is restarted as well
     *       causing unwanted disturbances on that channel.
     */
    g_awg_reg->state_machine_conf = 0x110011;

    fpga_awg_exit();
}
