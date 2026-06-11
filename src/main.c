#include "pebble/pebble.h"
#include "pebble/pebble_bigwig.h"
#include "pebble/pebble_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    pebble_config_t config;
    const char *input_path;
    const char *output_path;
    const char *chrom_sizes_path;
    const char *chrom_filter;
    pebble_input_format_t format;
    pebble_output_format_t output_format;
    int demo_mode;
    int format_set;
    int output_format_set;
    int want_bigwig;
    int no_smooth;
} cli_options_t;

static int path_ends_with(const char *path, const char *suffix)
{
    size_t path_len;
    size_t suffix_len;

    if (path == NULL || suffix == NULL) {
        return 0;
    }

    path_len = strlen(path);
    suffix_len = strlen(suffix);
    if (path_len < suffix_len) {
        return 0;
    }
    return strcmp(path + path_len - suffix_len, suffix) == 0;
}

static char *replace_path_suffix(const char *path, const char *old_suffix, const char *new_suffix)
{
    size_t path_len = strlen(path);
    size_t old_len = strlen(old_suffix);
    size_t new_len = strlen(new_suffix);
    char *out;

    if (path_len < old_len) {
        return NULL;
    }
    out = (char *)malloc(path_len - old_len + new_len + 1U);
    if (out == NULL) {
        return NULL;
    }
    memcpy(out, path, path_len - old_len);
    memcpy(out + path_len - old_len, new_suffix, new_len + 1U);
    return out;
}

static char *derive_bedgraph_path(const char *output_path)
{
    if (output_path == NULL) {
        return NULL;
    }
    if (path_ends_with(output_path, ".bw")) {
        return replace_path_suffix(output_path, ".bw", ".bedgraph");
    }
    if (path_ends_with(output_path, ".bigwig")) {
        return replace_path_suffix(output_path, ".bigwig", ".bedgraph");
    }
    return strdup(output_path);
}

static char *derive_smoothed_bedgraph_path(const char *output_path)
{
    if (output_path == NULL) {
        return NULL;
    }
    if (path_ends_with(output_path, ".smoothed.bedgraph")) {
        return strdup(output_path);
    }
    if (path_ends_with(output_path, ".normalised.bedgraph")) {
        return replace_path_suffix(output_path, ".normalised.bedgraph", ".smoothed.bedgraph");
    }
    if (path_ends_with(output_path, ".bedgraph")) {
        return replace_path_suffix(output_path, ".bedgraph", ".smoothed.bedgraph");
    }
    if (path_ends_with(output_path, ".bw") || path_ends_with(output_path, ".bigwig")) {
        char *bedgraph = derive_bedgraph_path(output_path);
        char *smoothed;

        if (bedgraph == NULL) {
            return NULL;
        }
        smoothed = replace_path_suffix(bedgraph, ".bedgraph", ".smoothed.bedgraph");
        free(bedgraph);
        return smoothed;
    }
    return NULL;
}

static char *derive_normalised_bedgraph_path(const char *output_path)
{
    if (output_path == NULL) {
        return NULL;
    }
    if (path_ends_with(output_path, ".normalised.bedgraph")) {
        return strdup(output_path);
    }
    if (path_ends_with(output_path, ".smoothed.bedgraph")) {
        return replace_path_suffix(output_path, ".smoothed.bedgraph", ".normalised.bedgraph");
    }
    return derive_bedgraph_path(output_path);
}

static char *derive_bigwig_path(const char *output_path)
{
    char *out;

    if (output_path == NULL) {
        return NULL;
    }
    if (path_ends_with(output_path, ".bw") || path_ends_with(output_path, ".bigwig")) {
        return strdup(output_path);
    }
    if (path_ends_with(output_path, ".bedgraph")) {
        return replace_path_suffix(output_path, ".bedgraph", ".bw");
    }

    out = (char *)malloc(strlen(output_path) + 4U);
    if (out == NULL) {
        return NULL;
    }
    sprintf(out, "%s.bw", output_path);
    return out;
}

static void print_usage(const char *program)
{
    fprintf(
        stderr,
        "Usage: %s [--demo] [-i input] [-o output] [-c chrom] [-f bed|bedgraph]\n"
        "              [-W window] [-S step] [--trim-low frac] [--trim-high frac]\n"
        "              [--sizes chrom.sizes] [--format-out bedgraph|bigwig]\n"
        "              [--no-smooth]\n"
        "\n"
        "  --demo            Run the built-in mock scaffold demo\n"
        "  -i input           Input BED or BedGraph file\n"
        "  -o output          Output BedGraph file (default: stdout)\n"
        "  --sizes file      Chromosome sizes (name<TAB>length); also writes BigWig\n"
        "  --no-smooth       Convert coverage to BigWig without smoothing (needs --sizes, -o)\n"
        "  --format-out fmt  With --sizes: force BigWig basename (default: from -o)\n"
        "                    BigWig requires a build with -DPEBBLE_BIGWIG=ON\n"
        "  -c chrom           Process only the named chromosome or contig\n"
        "  -f bed|bedgraph    Input format (default: inferred from file extension)\n"
        "  -W window          Sliding window size (default: 1000)\n"
        "  -S step            Step size (default: 100)\n"
        "  --trim-low frac    Trim fraction from low end (default: 0.40)\n"
        "  --trim-high frac   Trim fraction from high end (default: 0.40)\n"
        "\n"
        "BedGraph is always written. With --sizes, a .bw file is also written\n"
        "  (-o out.bedgraph -> out.bedgraph + out.bw; -o out.bw -> out.bw + out.bedgraph).\n"
        "\n"
        "Examples:\n"
        "  %s -i coverage.bedgraph -o smoothed.bedgraph\n"
        "  %s -i coverage.bed -o smoothed.bedgraph --sizes chrom.sizes\n"
        "  %s -i coverage.bedgraph -o out.bw --sizes chrom.sizes --no-smooth\n",
        program,
        program,
        program,
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

static int parse_output_format(const char *text, pebble_output_format_t *format)
{
    if (strcmp(text, "bedgraph") == 0) {
        *format = PEBBLE_OUTPUT_BEDGRAPH;
        return 0;
    }
    if (strcmp(text, "bigwig") == 0 || strcmp(text, "bw") == 0) {
        *format = PEBBLE_OUTPUT_BIGWIG;
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
    options->chrom_sizes_path = NULL;
    options->chrom_filter = NULL;
    options->format = PEBBLE_FORMAT_BEDGRAPH;
    options->output_format = PEBBLE_OUTPUT_BEDGRAPH;
    options->demo_mode = 0;
    options->format_set = 0;
    options->output_format_set = 0;
    options->want_bigwig = 0;
    options->no_smooth = 0;

    for (int i = 1; i < argc; i++) {
        const char *arg = argv[i];

        if (strcmp(arg, "--demo") == 0) {
            options->demo_mode = 1;
            continue;
        }
        if (strcmp(arg, "--no-smooth") == 0) {
            options->no_smooth = 1;
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
        if (strcmp(arg, "--sizes") == 0 || strcmp(arg, "-sizes") == 0) {
            if (i + 1 >= argc) {
                return -1;
            }
            options->chrom_sizes_path = argv[++i];
            continue;
        }
        if (strcmp(arg, "--format-out") == 0) {
            if (i + 1 >= argc) {
                return -1;
            }
            if (parse_output_format(argv[++i], &options->output_format) != 0) {
                return -1;
            }
            options->output_format_set = 1;
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
    if (options->output_path != NULL && !options->output_format_set) {
        options->output_format = pebble_infer_output_format(options->output_path);
    }

    options->want_bigwig = options->chrom_sizes_path != NULL
        || options->output_format == PEBBLE_OUTPUT_BIGWIG;

    if (options->no_smooth) {
        options->want_bigwig = 1;
    }

    if (options->want_bigwig && options->output_path == NULL) {
        return -1;
    }
    if (options->want_bigwig && options->chrom_sizes_path == NULL) {
        return -1;
    }
    if (options->no_smooth && options->demo_mode) {
        return -1;
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

    {
        double weighted_sum = 0.0;
        size_t total_bases = 0;
        double genome_average = 0.0;
        double normalise_factor = 1.0;

        pebble_smoothed_genome_average_add(0, config, smoothed_track, out_len, &weighted_sum, &total_bases);
        if (total_bases > 0U) {
            genome_average = weighted_sum / (double)total_bases;
            normalise_factor = pebble_coverage_normalise_factor(genome_average);
            fprintf(
                stderr,
                "pebble: genome average=%.2f, normalise factor=%.4f (avg/2), target mean=2.0\n",
                genome_average,
                normalise_factor
            );
        }

        printf("Interval_Start\tInterval_End\tSmoothed\tNormalised\n");
        for (size_t idx = 0; idx < out_len; idx++) {
            int start = 0;
            int end = 0;
            double smoothed_value = smoothed_track[idx];
            double normalised_value = smoothed_value / normalise_factor;

            pebble_output_interval(idx, 0, config, &start, &end);
            printf(
                "%d\t%d\t%d\t%d\n",
                start,
                end,
                pebble_round_coverage(smoothed_value),
                pebble_round_coverage(normalised_value)
            );
        }
    }

    free(smoothed_track);
    return EXIT_SUCCESS;
}

typedef struct {
    const pebble_coverage_t *entry;
    double *values;
    size_t value_count;
} smoothed_contig_t;

static void free_smoothed_contigs(smoothed_contig_t *contigs, size_t count)
{
    size_t i;

    if (contigs == NULL) {
        return;
    }

    for (i = 0; i < count; i++) {
        free(contigs[i].values);
    }
    free(contigs);
}

static int process_batch(
    const pebble_coverage_batch_t *batch,
    const cli_options_t *options,
    FILE *smoothed_output,
    FILE *normalised_output,
    const pebble_chrom_sizes_t *chrom_sizes,
    int validate_chrom_sizes,
    size_t *processed_out,
    size_t *skipped_out)
{
    smoothed_contig_t *contigs = NULL;
    size_t contig_count = 0;
    size_t processed = 0;
    size_t skipped = 0;
    double weighted_sum = 0.0;
    size_t total_bases = 0;
    double genome_average = 0.0;
    double normalise_factor = 1.0;

    contigs = (smoothed_contig_t *)calloc(batch->count, sizeof(smoothed_contig_t));
    if (contigs == NULL) {
        fprintf(stderr, "pebble: out of memory\n");
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < batch->count; i++) {
        const pebble_coverage_t *entry = &batch->items[i];
        size_t num_steps = pebble_output_count(entry->length, &options->config);
        size_t out_len = 0;
        double *smoothed_track = NULL;
        pebble_status_t status;

        if (num_steps == 0U) {
            fprintf(
                stderr,
                "pebble: skipping %s (length %zu shorter than window %d)\n",
                entry->chrom,
                entry->length,
                options->config.window_size
            );
            skipped++;
            continue;
        }

        if (validate_chrom_sizes && !pebble_chrom_sizes_contains(chrom_sizes, entry->chrom)) {
            fprintf(
                stderr,
                "pebble: %s not found in chromosome sizes file\n",
                entry->chrom
            );
            free_smoothed_contigs(contigs, contig_count);
            if (processed_out != NULL) {
                *processed_out = processed;
            }
            if (skipped_out != NULL) {
                *skipped_out = skipped;
            }
            return EXIT_FAILURE;
        }

        smoothed_track = (double *)malloc(num_steps * sizeof(double));
        if (smoothed_track == NULL) {
            fprintf(stderr, "pebble: out of memory\n");
            free_smoothed_contigs(contigs, contig_count);
            if (processed_out != NULL) {
                *processed_out = processed;
            }
            if (skipped_out != NULL) {
                *skipped_out = skipped;
            }
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
            free_smoothed_contigs(contigs, contig_count);
            if (processed_out != NULL) {
                *processed_out = processed;
            }
            if (skipped_out != NULL) {
                *skipped_out = skipped;
            }
            return EXIT_FAILURE;
        }

        pebble_smoothed_genome_average_add(
            entry->start_offset,
            &options->config,
            smoothed_track,
            out_len,
            &weighted_sum,
            &total_bases
        );

        contigs[contig_count].entry = entry;
        contigs[contig_count].values = smoothed_track;
        contigs[contig_count].value_count = out_len;
        contig_count++;
        processed++;
    }

    if (total_bases > 0U) {
        genome_average = weighted_sum / (double)total_bases;
        normalise_factor = pebble_coverage_normalise_factor(genome_average);
        fprintf(
            stderr,
            "pebble: genome average=%.2f, normalise factor=%.4f (avg/2), target mean=2.0\n",
            genome_average,
            normalise_factor
        );
    }

    for (size_t i = 0; i < contig_count; i++) {
        const pebble_coverage_t *entry = contigs[i].entry;
        pebble_io_status_t io_status;

        if (smoothed_output != NULL) {
            io_status = pebble_write_bedgraph(
                smoothed_output,
                entry->chrom,
                entry->start_offset,
                &options->config,
                contigs[i].values,
                contigs[i].value_count
            );
            if (io_status != PEBBLE_IO_OK) {
                free_smoothed_contigs(contigs, contig_count);
                if (processed_out != NULL) {
                    *processed_out = processed;
                }
                if (skipped_out != NULL) {
                    *skipped_out = skipped;
                }
                fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
                return EXIT_FAILURE;
            }
        }

        pebble_normalise_smoothed_values(
            contigs[i].values,
            contigs[i].value_count,
            normalise_factor
        );

        if (normalised_output != NULL) {
            io_status = pebble_write_bedgraph(
                normalised_output,
                entry->chrom,
                entry->start_offset,
                &options->config,
                contigs[i].values,
                contigs[i].value_count
            );
            if (io_status != PEBBLE_IO_OK) {
                free_smoothed_contigs(contigs, contig_count);
                if (processed_out != NULL) {
                    *processed_out = processed;
                }
                if (skipped_out != NULL) {
                    *skipped_out = skipped;
                }
                fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
                return EXIT_FAILURE;
            }
        }
    }

    free_smoothed_contigs(contigs, contig_count);

    if (processed_out != NULL) {
        *processed_out = processed;
    }
    if (skipped_out != NULL) {
        *skipped_out = skipped;
    }
    return EXIT_SUCCESS;
}

static int write_bigwig_from_bedgraph_file(
    const char *bedgraph_path,
    const char *bigwig_path,
    const pebble_chrom_sizes_t *chrom_sizes)
{
    pebble_bedgraph_batch_t records = {0};
    pebble_bigwig_writer_t *writer = NULL;
    pebble_io_status_t io_status;

    io_status = pebble_read_bedgraph_records(bedgraph_path, NULL, &records);
    if (io_status != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: failed to read BedGraph for BigWig (%s)\n", pebble_io_status_string(io_status));
        return EXIT_FAILURE;
    }

    pebble_bedgraph_batch_sort(&records);

    writer = pebble_bigwig_writer_create(bigwig_path, chrom_sizes);
    if (writer == NULL) {
        fprintf(stderr, "pebble: failed to create BigWig output file\n");
        pebble_bedgraph_batch_free(&records);
        return EXIT_FAILURE;
    }

    fprintf(
        stderr,
        "pebble: writing BigWig from BedGraph (%zu interval(s), same as bedGraphToBigWig)\n",
        records.count
    );

    io_status = pebble_bigwig_write_records(writer, records.items, records.count);
    if (io_status != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
        pebble_bigwig_writer_close(writer);
        pebble_bedgraph_batch_free(&records);
        return EXIT_FAILURE;
    }

    if (pebble_bigwig_writer_close(writer) != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: failed to finalize BigWig output\n");
        pebble_bedgraph_batch_free(&records);
        return EXIT_FAILURE;
    }

    pebble_bedgraph_batch_free(&records);
    return EXIT_SUCCESS;
}

static int run_file_mode(const cli_options_t *options)
{
    pebble_coverage_batch_t batch = {0};
    pebble_chrom_sizes_t chrom_sizes = {0};
    pebble_io_status_t io_status;
    FILE *smoothed_output = NULL;
    FILE *normalised_output = stdout;
    char *smoothed_path = NULL;
    char *normalised_path = NULL;
    char *bigwig_path = NULL;
    int exit_code = EXIT_SUCCESS;
    size_t processed = 0;
    size_t skipped = 0;

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

    fprintf(
        stderr,
        "pebble: loaded %zu contig(s) from %s\n",
        batch.count,
        options->input_path
    );

    pebble_coverage_batch_sort(&batch);

    if (options->output_path != NULL) {
        normalised_path = derive_normalised_bedgraph_path(options->output_path);
        smoothed_path = derive_smoothed_bedgraph_path(options->output_path);
        if (normalised_path == NULL || smoothed_path == NULL) {
            fprintf(stderr, "pebble: out of memory\n");
            free(normalised_path);
            free(smoothed_path);
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }

        normalised_output = fopen(normalised_path, "w");
        if (normalised_output == NULL) {
            fprintf(stderr, "pebble: failed to open normalised BedGraph output file\n");
            free(normalised_path);
            free(smoothed_path);
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }

        smoothed_output = fopen(smoothed_path, "w");
        if (smoothed_output == NULL) {
            fprintf(stderr, "pebble: failed to open smoothed BedGraph output file\n");
            fclose(normalised_output);
            free(normalised_path);
            free(smoothed_path);
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }
    }

    if (options->want_bigwig) {
        bigwig_path = derive_bigwig_path(options->output_path);
        if (bigwig_path == NULL) {
            fprintf(stderr, "pebble: out of memory\n");
            if (smoothed_path != NULL) {
                fclose(smoothed_output);
                fclose(normalised_output);
                free(smoothed_path);
                free(normalised_path);
            }
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }

        io_status = pebble_read_chrom_sizes(options->chrom_sizes_path, &chrom_sizes);
        if (io_status != PEBBLE_IO_OK) {
            fprintf(
                stderr,
                "pebble: failed to read chromosome sizes (%s)\n",
                pebble_io_status_string(io_status)
            );
            if (smoothed_path != NULL) {
                fclose(smoothed_output);
                fclose(normalised_output);
                free(smoothed_path);
                free(normalised_path);
            }
            free(bigwig_path);
            pebble_coverage_batch_free(&batch);
            return EXIT_FAILURE;
        }

        fprintf(
            stderr,
            "pebble: writing smoothed BedGraph to %s, normalised BedGraph to %s; BigWig from normalised -> %s\n",
            smoothed_path,
            normalised_path,
            bigwig_path
        );
    } else if (normalised_path != NULL) {
        fprintf(
            stderr,
            "pebble: writing smoothed BedGraph to %s, normalised BedGraph to %s\n",
            smoothed_path,
            normalised_path
        );
    } else {
        fprintf(stderr, "pebble: writing normalised BedGraph to stdout (smoothed not written)\n");
    }

    exit_code = process_batch(
        &batch,
        options,
        smoothed_output,
        normalised_output,
        &chrom_sizes,
        options->want_bigwig,
        &processed,
        &skipped
    );

    if (smoothed_path != NULL) {
        if (fclose(smoothed_output) != 0 && exit_code == EXIT_SUCCESS) {
            fprintf(stderr, "pebble: failed to close smoothed BedGraph output file\n");
            exit_code = EXIT_FAILURE;
        }
        if (fclose(normalised_output) != 0 && exit_code == EXIT_SUCCESS) {
            fprintf(stderr, "pebble: failed to close normalised BedGraph output file\n");
            exit_code = EXIT_FAILURE;
        }
    }

    if (exit_code == EXIT_SUCCESS && options->want_bigwig && normalised_path != NULL) {
        exit_code = write_bigwig_from_bedgraph_file(normalised_path, bigwig_path, &chrom_sizes);
    }

    fprintf(
        stderr,
        "pebble: done — processed %zu contig(s), skipped %zu (window %d)\n",
        processed,
        skipped,
        options->config.window_size
    );

    if (exit_code == EXIT_SUCCESS && processed == 0U) {
        fprintf(
            stderr,
            "pebble: no output written (all contigs shorter than window, or none loaded)\n"
        );
        exit_code = EXIT_FAILURE;
    }

    free(smoothed_path);
    free(normalised_path);
    free(bigwig_path);
    pebble_chrom_sizes_free(&chrom_sizes);
    pebble_coverage_batch_free(&batch);
    return exit_code;
}

static int run_convert_mode(const cli_options_t *options)
{
    pebble_bedgraph_batch_t batch = {0};
    pebble_chrom_sizes_t chrom_sizes = {0};
    pebble_bigwig_writer_t *writer = NULL;
    char *bigwig_path = NULL;
    pebble_io_status_t io_status;
    int exit_code = EXIT_SUCCESS;
    size_t i;

    io_status = pebble_read_interval_records(
        options->input_path,
        options->chrom_filter,
        options->format == PEBBLE_FORMAT_BED ? 1 : 0,
        &batch
    );
    if (io_status != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
        return EXIT_FAILURE;
    }

    fprintf(
        stderr,
        "pebble: loaded %zu BedGraph interval(s) from %s\n",
        batch.count,
        options->input_path
    );

    pebble_bedgraph_batch_sort(&batch);

    bigwig_path = derive_bigwig_path(options->output_path);
    if (bigwig_path == NULL) {
        fprintf(stderr, "pebble: out of memory\n");
        pebble_bedgraph_batch_free(&batch);
        return EXIT_FAILURE;
    }

    io_status = pebble_read_chrom_sizes(options->chrom_sizes_path, &chrom_sizes);
    if (io_status != PEBBLE_IO_OK) {
        fprintf(
            stderr,
            "pebble: failed to read chromosome sizes (%s)\n",
            pebble_io_status_string(io_status)
        );
        free(bigwig_path);
        pebble_bedgraph_batch_free(&batch);
        return EXIT_FAILURE;
    }

    for (i = 0; i < batch.count; i++) {
        if (!pebble_chrom_sizes_contains(&chrom_sizes, batch.items[i].chrom)) {
            fprintf(
                stderr,
                "pebble: %s not found in chromosome sizes file\n",
                batch.items[i].chrom
            );
            free(bigwig_path);
            pebble_chrom_sizes_free(&chrom_sizes);
            pebble_bedgraph_batch_free(&batch);
            return EXIT_FAILURE;
        }
    }

    writer = pebble_bigwig_writer_create(bigwig_path, &chrom_sizes);
    if (writer == NULL) {
        fprintf(stderr, "pebble: failed to create BigWig output file\n");
        free(bigwig_path);
        pebble_chrom_sizes_free(&chrom_sizes);
        pebble_bedgraph_batch_free(&batch);
        return EXIT_FAILURE;
    }

    fprintf(
        stderr,
        "pebble: converting to BigWig (no smoothing, bedGraphToBigWig-style) -> %s\n",
        bigwig_path
    );

    io_status = pebble_bigwig_write_records(writer, batch.items, batch.count);
    if (io_status != PEBBLE_IO_OK) {
        fprintf(stderr, "pebble: %s\n", pebble_io_status_string(io_status));
        exit_code = EXIT_FAILURE;
    }

    if (pebble_bigwig_writer_close(writer) != PEBBLE_IO_OK && exit_code == EXIT_SUCCESS) {
        fprintf(stderr, "pebble: failed to finalize BigWig output\n");
        exit_code = EXIT_FAILURE;
    }

    if (exit_code == EXIT_SUCCESS) {
        fprintf(stderr, "pebble: done — wrote %zu interval(s) to BigWig\n", batch.count);
    }

    free(bigwig_path);
    pebble_chrom_sizes_free(&chrom_sizes);
    pebble_bedgraph_batch_free(&batch);
    return exit_code;
}

static int bigwig_unavailable(const cli_options_t *options)
{
#ifndef PEBBLE_BIGWIG
    if (options->want_bigwig || options->chrom_sizes_path != NULL || options->no_smooth) {
        fprintf(
            stderr,
            "pebble: BigWig output is not enabled in this build "
            "(rebuild with -DPEBBLE_BIGWIG=ON or make PEBBLE_BIGWIG=1)\n"
        );
        return 1;
    }
#else
    (void)options;
#endif
    return 0;
}

int main(int argc, char **argv)
{
    cli_options_t options;

    if (parse_options(argc, argv, &options) != 0) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (bigwig_unavailable(&options)) {
        return EXIT_FAILURE;
    }

    if (options.demo_mode) {
        return run_demo(&options.config);
    }

    if (options.no_smooth) {
        return run_convert_mode(&options);
    }

    return run_file_mode(&options);
}
