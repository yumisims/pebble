#include "pebble/pebble.h"
#include "pebble/pebble_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    pebble_config_t config;
    const char *input_path;
    const char *output_path;
    const char *chrom_filter;
    pebble_input_format_t format;
    int demo_mode;
    int format_set;
} cli_options_t;

static void print_usage(const char *program)
{
    fprintf(
        stderr,
        "Usage: %s [--demo] [-i input] [-o output] [-c chrom] [-f bed|bedgraph]\n"
        "              [-W window] [-S step] [--trim-low frac] [--trim-high frac]\n"
        "\n"
        "  --demo            Run the built-in mock scaffold demo\n"
        "  -i input           Input BED or BedGraph file\n"
        "  -o output          Output BedGraph file (default: stdout)\n"
        "  -c chrom           Process only the named chromosome or contig\n"
        "  -f bed|bedgraph    Input format (default: inferred from file extension)\n"
        "  -W window          Sliding window size (default: 1000)\n"
        "  -S step            Step size (default: 100)\n"
        "  --trim-low frac    Trim fraction from low end (default: 0.40)\n"
        "  --trim-high frac   Trim fraction from high end (default: 0.40)\n",
        program
    );
}

static pebble_input_format_t infer_format(const char *path)
{
    size_t len = strlen(path);

    if (len >= 9 && strcmp(path + len - 9, ".bedgraph") == 0) {
        return PEBBLE_FORMAT_BEDGRAPH;
    }
    if (len >= 4 && strcmp(path + len - 4, ".bed") == 0) {
        return PEBBLE_FORMAT_BED;
    }
    return PEBBLE_FORMAT_BEDGRAPH;
}

static int parse_format(const char *text, pebble_input_format_t *format)
{
    if (strcmp(text, "bedgraph") == 0) {
        *format = PEBBLE_FORMAT_BEDGRAPH;
        return 0;
    }
    if (strcmp(text, "bed") == 0) {
        *format = PEBBLE_FORMAT_BED;
        return 0;
    }
    return -1;
}

static int option_requires_value(char opt)
{
    switch (opt) {
    case 'i':
    case 'o':
    case 'c':
    case 'f':
    case 'W':
    case 'S':
        return 1;
    default:
        return 0;
    }
}

static int parse_options(int argc, char **argv, cli_options_t *options)
{
    options->config.window_size = 1000;
    options->config.step_size = 100;
    options->config.trim_low = 0.40;
    options->config.trim_high = 0.40;
    options->input_path = NULL;
    options->output_path = NULL;
    options->chrom_filter = NULL;
    options->format = PEBBLE_FORMAT_BEDGRAPH;
    options->demo_mode = 0;
    options->format_set = 0;

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];

        if (strcmp(arg, "--demo") == 0) {
            options->demo_mode = 1;
            continue;
        }
        if (strcmp(arg, "--trim-low") == 0) {
            if (i + 1 >= argc) {
                return -1;
            }
            options->config.trim_low = atof(argv[++i]);
            continue;
        }
        if (strcmp(arg, "--trim-high") == 0) {
            if (i + 1 >= argc) {
                return -1;
            }
            options->config.trim_high = atof(argv[++i]);
            continue;
        }
        if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
            return -1;
        }
        if (arg[0] != '-' || arg[1] == '\0') {
            return -1;
        }

        for (int j = 1; arg[j] != '\0'; j++) {
            char opt = arg[j];
            const char *val = NULL;

            if (option_requires_value(opt)) {
                if (arg[j + 1] != '\0') {
                    val = arg + j + 1;
                    j = (int)strlen(arg) - 1;
                } else {
                    if (i + 1 >= argc) {
                        return -1;
                    }
                    val = argv[++i];
                }
            }

            switch (opt) {
            case 'i':
                options->input_path = val;
                break;
            case 'o':
                options->output_path = val;
                break;
            case 'c':
                options->chrom_filter = val;
                break;
            case 'f':
                if (val == NULL || parse_format(val, &options->format) != 0) {
                    return -1;
                }
                options->format_set = 1;
                break;
            case 'W':
                if (val == NULL) {
                    return -1;
                }
                options->config.window_size = atoi(val);
                break;
            case 'S':
                if (val == NULL) {
                    return -1;
                }
                options->config.step_size = atoi(val);
                break;
            case 'h':
                return -1;
            default:
                return -1;
            }
        }
    }

    if (options->input_path == NULL && options->demo_mode == 0) {
        options->demo_mode = 1;
    }
    if (options->input_path != NULL && !options->format_set) {
        options->format = infer_format(options->input_path);
    }
    return 0;
}

static void print_status(pebble_status_t status)
{
    switch (status) {
    case PEBBLE_OK:
        break;
    case PEBBLE_ERR_INVALID_ARG:
        fprintf(stderr, "pebble: invalid argument\n");
        break;
    case PEBBLE_ERR_TOO_SHORT:
        fprintf(stderr, "pebble: coverage shorter than window size\n");
        break;
    case PEBBLE_ERR_BUFFER_TOO_SMALL:
        fprintf(stderr, "pebble: output buffer too small\n");
        break;
    default:
        fprintf(stderr, "pebble: unknown error (%d)\n", status);
        break;
    }
}

static int run_demo(const pebble_config_t *config)
{
    int16_t mock_coverage[1500];
    size_t out_len = 0;
    size_t num_steps = 0;
    double *smoothed_track = NULL;
    pebble_status_t status;

    for (size_t i = 0; i < 1500; i++) {
        mock_coverage[i] = 30;
    }
    for (size_t i = 200; i < 250; i++) {
        mock_coverage[i] = 0;
    }
    for (size_t i = 700; i < 850; i++) {
        mock_coverage[i] = 999;
    }

    num_steps = pebble_output_count(1500, config);
    smoothed_track = (double *)malloc(num_steps * sizeof(double));
    if (smoothed_track == NULL) {
        fprintf(stderr, "pebble: out of memory\n");
        return EXIT_FAILURE;
    }

    status = pebble_process(mock_coverage, 1500, config, smoothed_track, num_steps, &out_len);
    if (status != PEBBLE_OK) {
        print_status(status);
        free(smoothed_track);
        return EXIT_FAILURE;
    }

    printf("Interval_Start\tInterval_End\tSmoothed_Value\n");
    for (size_t idx = 0; idx < out_len; idx++) {
        int start = 0;
        int end = 0;

        pebble_output_interval(idx, 0, config, &start, &end);
        printf("%d\t%d\t%d\n", start, end, pebble_round_coverage(smoothed_track[idx]));
    }

    free(smoothed_track);
    return EXIT_SUCCESS;
}

static int process_batch(
    const pebble_coverage_batch_t *batch,
    const cli_options_t *options,
    FILE *output)
{
    for (size_t i = 0; i < batch->count; i++) {
        const pebble_coverage_t *entry = &batch->items[i];
        size_t num_steps = pebble_output_count(entry->length, &options->config);
        size_t out_len = 0;
        double *smoothed_track = NULL;
        pebble_status_t status;
        pebble_io_status_t io_status;

        if (num_steps == 0U) {
            fprintf(
                stderr,
                "pebble: skipping %s (length %zu shorter than window %d)\n",
                entry->chrom,
                entry->length,
                options->config.window_size
            );
            continue;
        }

        smoothed_track = (double *)malloc(num_steps * sizeof(double));
        if (smoothed_track == NULL) {
            fprintf(stderr, "pebble: out of memory\n");
            return EXIT_FAILURE;
        }

        status = pebble_process(
            entry->coverage,
            entry->length,
            &options->config,
            smoothed_track,
            num_steps,
            &out_len
        );
        if (status != PEBBLE_OK) {
            print_status(status);
            free(smoothed_track);
            return EXIT_FAILURE;
        }

        io_status = pebble_write_bedgraph(
            output,
            entry->chrom,
            entry->start_offset,
            &options->config,
            smoothed_track,
            out_len
        );
        free(smoothed_track);
        if (io_status != PEBBLE_IO_OK) {
            fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static int run_file_mode(const cli_options_t *options)
{
    pebble_coverage_batch_t batch = {0};
    pebble_io_status_t io_status;
    FILE *output = stdout;
    int exit_code = EXIT_SUCCESS;

    io_status = pebble_read_coverage(
        options->input_path,
        options->format,
        options->chrom_filter,
        &batch
    );
    if (io_status != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
        return EXIT_FAILURE;
    }

    if (options->output_path != NULL) {
        output = fopen(options->output_path, "w");
        if (output == NULL) {
            fprintf(stderr, "pebble: failed to open output file\n");
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }
    }

    exit_code = process_batch(&batch, options, output);

    if (options->output_path != NULL) {
        if (fclose(output) != 0 && exit_code == EXIT_SUCCESS) {
            fprintf(stderr, "pebble: failed to close output file\n");
            exit_code = EXIT_FAILURE;
        }
    }

    pebble_coverage_batch_free(&batch);
    return exit_code;
}

int main(int argc, char **argv)
{
    cli_options_t options;

    if (parse_options(argc, argv, &options) != 0) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (options.demo_mode) {
        return run_demo(&options.config);
    }

    return run_file_mode(&options);
}
