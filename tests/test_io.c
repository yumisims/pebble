#include "pebble/pebble.h"
#include "pebble/pebble_bigwig.h"
#include "pebble/pebble_io.h"

#ifdef PEBBLE_BIGWIG
#include "bigWig.h"
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_read_example_cov_bed(void)
{
    pebble_coverage_batch_t batch = {0};
    pebble_io_status_t status = pebble_read_bed(
        "examples/example_cov.bed",
        "scaffold_1",
        &batch
    );

    assert(status == PEBBLE_IO_OK);
    assert(batch.count == 1U);
    assert(strcmp(batch.items[0].chrom, "scaffold_1") == 0);
    assert(batch.items[0].start_offset == 0);
    assert(batch.items[0].length == 1527318U);
    assert(batch.items[0].coverage[0] == 4);
    assert(batch.items[0].coverage[559] == 5);

    pebble_coverage_batch_free(&batch);
}

static void test_batch_sort(void)
{
    pebble_coverage_batch_t batch = {0};

    batch.count = 2U;
    batch.items = (pebble_coverage_t *)calloc(2U, sizeof(pebble_coverage_t));
    assert(batch.items != NULL);
    batch.items[0].chrom = strdup("chrZ");
    batch.items[0].start_offset = 100;
    batch.items[1].chrom = strdup("chrA");
    batch.items[1].start_offset = 0;

    pebble_coverage_batch_sort(&batch);
    assert(strcmp(batch.items[0].chrom, "chrA") == 0);
    assert(strcmp(batch.items[1].chrom, "chrZ") == 0);

    pebble_coverage_batch_free(&batch);
}

static void test_gap_fill_at_start(void)
{
    pebble_config_t config = {
        .window_size = 100,
        .step_size = 50,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    double *out = NULL;
    size_t out_len = 0;
    size_t num_steps = 0;

    assert(pebble_read_bed("examples/gap_start.bed", "scaffold_1", &batch) == PEBBLE_IO_OK);
    assert(batch.count == 1U);
    assert(batch.items[0].start_offset == 0);
    assert(batch.items[0].coverage[0] == 5);
    assert(batch.items[0].coverage[999] == 5);
    assert(batch.items[0].coverage[1000] == 5);
    assert(batch.items[0].coverage[1100] == 5);

    num_steps = pebble_output_count(batch.items[0].length, &config);
    out = (double *)malloc(num_steps * sizeof(double));
    assert(out != NULL);
    assert(
        pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, num_steps, &out_len)
        == PEBBLE_OK
    );
    assert(out_len >= 1U);
    assert(out[0] >= 4.5);

    free(out);
    pebble_coverage_batch_free(&batch);
}

static void test_write_bedgraph_roundtrip(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    double *out = NULL;
    size_t out_len = 0;
    size_t num_steps = 0;
    FILE *tmp = tmpfile();
    char line[256];

    assert(pebble_read_bed("examples/example_cov.bed", "scaffold_1", &batch) == PEBBLE_IO_OK);
    num_steps = pebble_output_count(batch.items[0].length, &config);
    assert(num_steps > 0U);
    out = (double *)malloc(num_steps * sizeof(double));
    assert(out != NULL);
    assert(
        pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, num_steps, &out_len)
        == PEBBLE_OK
    );
    assert(out_len >= 1U);
    assert(pebble_write_bedgraph(tmp, "scaffold_1", 0, &config, out, out_len) == PEBBLE_IO_OK);

    rewind(tmp);
    assert(fgets(line, sizeof(line), tmp) != NULL);
    {
        char chrom[64];
        int start = -1;
        int end = -1;
        int value = -1;

        assert(sscanf(line, "%63s %d %d %d", chrom, &start, &end, &value) == 4);
        assert(strcmp(chrom, "scaffold_1") == 0);
        assert(start == 0);
        assert(end == 100);
        assert(value == 4);
    }

    free(out);
    fclose(tmp);
    pebble_coverage_batch_free(&batch);
}

static int bedgraph_fourth_column_is_integer(const char *line)
{
    const char *field = line;
    const char *next = NULL;
    size_t tab_count = 0;

    while (tab_count < 3U) {
        next = strchr(field, '\t');
        if (next == NULL) {
            return 0;
        }
        field = next + 1;
        tab_count++;
    }

    if (*field == '\0' || *field == '\n' || *field == '\r') {
        return 0;
    }

    while (*field != '\0' && *field != '\n' && *field != '\r') {
        if (*field == '.' || *field == 'e' || *field == 'E') {
            return 0;
        }
        field++;
    }

    return 1;
}

static void assert_bedgraph_integer_columns(FILE *handle)
{
    char line[256];

    rewind(handle);
    while (fgets(line, sizeof(line), handle) != NULL) {
        assert(bedgraph_fourth_column_is_integer(line));
    }
}

static void test_bedgraph_fourth_column_is_integer(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    double *smoothed = NULL;
    double *normalised = NULL;
    size_t out_len = 0;
    size_t num_steps = 0;
    double weighted_sum = 0.0;
    size_t total_bases = 0;
    double normalise_factor = 1.0;
    FILE *smoothed_tmp = tmpfile();
    FILE *normalised_tmp = tmpfile();

    assert(smoothed_tmp != NULL);
    assert(normalised_tmp != NULL);
    assert(pebble_read_bed("examples/example_cov.bed", "scaffold_1", &batch) == PEBBLE_IO_OK);
    num_steps = pebble_output_count(batch.items[0].length, &config);
    smoothed = (double *)malloc(num_steps * sizeof(double));
    normalised = (double *)malloc(num_steps * sizeof(double));
    assert(smoothed != NULL);
    assert(normalised != NULL);
    assert(
        pebble_process(batch.items[0].coverage, batch.items[0].length, &config, smoothed, num_steps, &out_len)
        == PEBBLE_OK
    );
    assert(out_len >= 1U);

    pebble_smoothed_genome_average_add(0, &config, smoothed, out_len, &weighted_sum, &total_bases);
    assert(total_bases > 0U);
    normalise_factor = pebble_coverage_normalise_factor(weighted_sum / (double)total_bases);

    memcpy(normalised, smoothed, out_len * sizeof(double));
    pebble_round_smoothed_values(smoothed, out_len);
    assert(pebble_write_bedgraph(smoothed_tmp, "scaffold_1", 0, &config, smoothed, out_len) == PEBBLE_IO_OK);
    assert_bedgraph_integer_columns(smoothed_tmp);

    pebble_normalise_smoothed_values(normalised, out_len, normalise_factor);
    assert(pebble_write_bedgraph(normalised_tmp, "scaffold_1", 0, &config, normalised, out_len) == PEBBLE_IO_OK);
    assert_bedgraph_integer_columns(normalised_tmp);

    free(smoothed);
    free(normalised);
    fclose(smoothed_tmp);
    fclose(normalised_tmp);
    pebble_coverage_batch_free(&batch);
}

static void test_coverage_normalisation(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    double *out = NULL;
    size_t out_len = 0;
    size_t num_steps = 0;
    double weighted_sum = 0.0;
    size_t total_bases = 0;
    double genome_average = 0.0;
    double normalise_factor = 0.0;
    double normalised_sum = 0.0;
    size_t normalised_bases = 0;

    assert(pebble_read_bed("examples/example_cov.bed", "scaffold_1", &batch) == PEBBLE_IO_OK);
    num_steps = pebble_output_count(batch.items[0].length, &config);
    out = (double *)malloc(num_steps * sizeof(double));
    assert(out != NULL);
    assert(
        pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, num_steps, &out_len)
        == PEBBLE_OK
    );

    pebble_smoothed_genome_average_add(0, &config, out, out_len, &weighted_sum, &total_bases);
    assert(total_bases > 0U);
    genome_average = weighted_sum / (double)total_bases;
    normalise_factor = pebble_coverage_normalise_factor(genome_average);
    assert(genome_average > 0.0);
    assert(normalise_factor > 0.0);

    pebble_normalise_smoothed_values(out, out_len, normalise_factor);
    pebble_smoothed_genome_average_add(0, &config, out, out_len, &normalised_sum, &normalised_bases);
    assert(normalised_bases > 0U);
    assert(fabs((normalised_sum / (double)normalised_bases) - 2.0) < 0.01);

    free(out);
    pebble_coverage_batch_free(&batch);
}

#ifdef PEBBLE_BIGWIG
static void test_write_bigwig_roundtrip(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    pebble_chrom_sizes_t sizes = {0};
    double *out = NULL;
    size_t out_len = 0;
    size_t num_steps = 0;
    const char *path = "build/test_example_cov.bw";
    pebble_bigwig_writer_t *writer = NULL;
    bigWigFile_t *reader = NULL;

    assert(pebble_read_chrom_sizes("examples/chrom.sizes", &sizes) == PEBBLE_IO_OK);
    assert(pebble_read_bed("examples/example_cov.bed", "scaffold_1", &batch) == PEBBLE_IO_OK);
    num_steps = pebble_output_count(batch.items[0].length, &config);
    out = (double *)malloc(num_steps * sizeof(double));
    assert(out != NULL);
    assert(
        pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, num_steps, &out_len)
        == PEBBLE_OK
    );
    assert(out_len >= 1U);
    assert(out[0] > 0.0);

    {
        FILE *bedgraph_tmp = fopen("build/test_example_cov.bedgraph", "w");
        pebble_bedgraph_batch_t records = {0};

        assert(bedgraph_tmp != NULL);
        assert(
            pebble_write_bedgraph(bedgraph_tmp, "scaffold_1", 0, &config, out, out_len) == PEBBLE_IO_OK
        );
        fclose(bedgraph_tmp);

        assert(pebble_read_bedgraph_records("build/test_example_cov.bedgraph", NULL, &records) == PEBBLE_IO_OK);
        assert(records.count == out_len);
        writer = pebble_bigwig_writer_create(path, &sizes);
        assert(writer != NULL);
        assert(pebble_bigwig_write_records(writer, records.items, records.count) == PEBBLE_IO_OK);
        assert(pebble_bigwig_writer_close(writer) == PEBBLE_IO_OK);
        pebble_bedgraph_batch_free(&records);
    }

    assert(bwInit(1U << 17) == 0);
    reader = bwOpen(path, NULL, "r");
    assert(reader != NULL);
    assert(bwIsBigWig(path, NULL) == 1);
    bwClose(reader);
    bwCleanup();

    free(out);
    pebble_chrom_sizes_free(&sizes);
    pebble_coverage_batch_free(&batch);
}

static void test_convert_bed_to_bigwig(void)
{
    pebble_bedgraph_batch_t batch = {0};
    pebble_chrom_sizes_t sizes = {0};
    const char *path = "build/test_convert.bw";
    pebble_bigwig_writer_t *writer = NULL;
    bigWigFile_t *reader = NULL;

    assert(pebble_read_chrom_sizes("examples/chrom.sizes", &sizes) == PEBBLE_IO_OK);
    assert(
        pebble_read_interval_records("examples/example_cov.bed", "scaffold_1", 1, &batch)
        == PEBBLE_IO_OK
    );
    assert(batch.count > 0U);
    assert(batch.items[0].start == 0);
    pebble_bedgraph_batch_sort(&batch);

    writer = pebble_bigwig_writer_create(path, &sizes);
    assert(writer != NULL);
    assert(pebble_bigwig_write_records(writer, batch.items, batch.count) == PEBBLE_IO_OK);
    assert(pebble_bigwig_writer_close(writer) == PEBBLE_IO_OK);

    assert(bwInit(1U << 17) == 0);
    reader = bwOpen(path, NULL, "r");
    assert(reader != NULL);
    assert(bwIsBigWig(path, NULL) == 1);
    bwClose(reader);
    bwCleanup();

    pebble_chrom_sizes_free(&sizes);
    pebble_bedgraph_batch_free(&batch);
}
#endif

int main(void)
{
    test_read_example_cov_bed();
    test_batch_sort();
    test_gap_fill_at_start();
    test_write_bedgraph_roundtrip();
    test_bedgraph_fourth_column_is_integer();
    test_coverage_normalisation();
#ifdef PEBBLE_BIGWIG
    test_write_bigwig_roundtrip();
    test_convert_bed_to_bigwig();
#endif

    puts("all io tests passed");
    return EXIT_SUCCESS;
}
