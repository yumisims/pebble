#include "pebble/pebble.h"
#include "pebble/pebble_io.h"

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
    assert(strncmp(line, "scaffold\t0\t100\t30.00", 20) == 0);

    fclose(tmp);
    pebble_coverage_batch_free(&batch);
}

int main(void)
{
    test_read_mock_bedgraph();
    test_read_mock_bed();
    test_read_genome2cov_bed();
    test_write_bedgraph_roundtrip();

    puts("all io tests passed");
    return EXIT_SUCCESS;
}
