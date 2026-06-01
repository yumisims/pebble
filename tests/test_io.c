#include "pebble/pebble.h"
#include "pebble/pebble_bigwig.h"
#include "pebble/pebble_io.h"

#ifdef PEBBLE_BIGWIG
#include "bigWig.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void test_read_mock_bedgraph(void)
{
    pebble_coverage_batch_t batch = {0};
    pebble_io_status_t status = pebble_read_bedgraph(
        "examples/mock_scaffold.bedgraph",
        "scaffold",
        &batch
    );

    assert(status == PEBBLE_IO_OK);
    assert(batch.count == 1U);
    assert(strcmp(batch.items[0].chrom, "scaffold") == 0);
    assert(batch.items[0].start_offset == 0);
    assert(batch.items[0].length == 1500U);
    assert(batch.items[0].coverage[0] == 30);
    assert(batch.items[0].coverage[225] == 0);
    assert(batch.items[0].coverage[800] == 999);

    pebble_coverage_batch_free(&batch);
}

static void test_read_mock_bed(void)
{
    pebble_coverage_batch_t batch = {0};
    pebble_io_status_t status = pebble_read_bed(
        "examples/mock_scaffold.bed",
        "scaffold",
        &batch
    );

    assert(status == PEBBLE_IO_OK);
    assert(batch.count == 1U);
    assert(batch.items[0].coverage[225] == 0);
    assert(batch.items[0].coverage[800] == 999);

    pebble_coverage_batch_free(&batch);
}

static void test_read_genome2cov_bed(void)
{
    pebble_coverage_batch_t batch = {0};
    pebble_io_status_t status = pebble_read_bed(
        "examples/genome2cov_style.bed",
        "HAP1_SCAFFOLD_1",
        &batch
    );

    assert(status == PEBBLE_IO_OK);
    assert(batch.count == 1U);
    assert(batch.items[0].coverage[0] == 2);
    assert(batch.items[0].coverage[93] == 42);

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

static void test_write_bedgraph_roundtrip(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };
    pebble_coverage_batch_t batch = {0};
    double out[6];
    size_t out_len = 0;
    FILE *tmp = tmpfile();
    char line[256];

    assert(pebble_read_bedgraph("examples/mock_scaffold.bedgraph", "scaffold", &batch) == PEBBLE_IO_OK);
    assert(pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, 6, &out_len) == PEBBLE_OK);
    assert(out_len == 6U);
    assert(pebble_write_bedgraph(tmp, "scaffold", 0, &config, out, out_len) == PEBBLE_IO_OK);

    rewind(tmp);
    assert(fgets(line, sizeof(line), tmp) != NULL);
    assert(strncmp(line, "scaffold\t0\t100\t30\n", 20) == 0);

    fclose(tmp);
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
    double out[6];
    size_t out_len = 0;
    const char *path = "build/test_mock_scaffold.bw";
    pebble_bigwig_writer_t *writer = NULL;
    bigWigFile_t *reader = NULL;

    assert(pebble_read_chrom_sizes("examples/mock_scaffold.sizes", &sizes) == PEBBLE_IO_OK);
    assert(pebble_read_bedgraph("examples/mock_scaffold.bedgraph", "scaffold", &batch) == PEBBLE_IO_OK);
    assert(pebble_process(batch.items[0].coverage, batch.items[0].length, &config, out, 6, &out_len) == PEBBLE_OK);
    assert(out_len == 6U);

    writer = pebble_bigwig_writer_create(path, &sizes);
    assert(writer != NULL);
    assert(
        pebble_bigwig_write_bedgraph(writer, "scaffold", 0, &config, out, out_len) == PEBBLE_IO_OK
    );
    assert(pebble_bigwig_writer_close(writer) == PEBBLE_IO_OK);

    assert(bwInit(1U << 17) == 0);
    reader = bwOpen(path, NULL, "r");
    assert(reader != NULL);
    assert(bwIsBigWig(path, NULL) == 1);
    bwClose(reader);
    bwCleanup();

    pebble_chrom_sizes_free(&sizes);
    pebble_coverage_batch_free(&batch);
}
#endif

#ifdef PEBBLE_BIGWIG
static void test_convert_bedgraph_to_bigwig(void)
{
    pebble_bedgraph_batch_t batch = {0};
    pebble_chrom_sizes_t sizes = {0};
    const char *path = "build/test_convert.bw";
    pebble_bigwig_writer_t *writer = NULL;
    bigWigFile_t *reader = NULL;

    assert(pebble_read_chrom_sizes("examples/mock_scaffold.sizes", &sizes) == PEBBLE_IO_OK);
    assert(
        pebble_read_bedgraph_records("examples/mock_scaffold.bedgraph", "scaffold", &batch)
        == PEBBLE_IO_OK
    );
    assert(batch.count == 3U);
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
    test_read_mock_bedgraph();
    test_read_mock_bed();
    test_read_genome2cov_bed();
    test_batch_sort();
    test_write_bedgraph_roundtrip();
#ifdef PEBBLE_BIGWIG
    test_write_bigwig_roundtrip();
    test_convert_bedgraph_to_bigwig();
#endif

    puts("all io tests passed");
    return EXIT_SUCCESS;
}
