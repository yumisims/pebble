#include "pebble/pebble_bigwig.h"

#include <stdlib.h>
#include <string.h>

void pebble_chrom_sizes_free(pebble_chrom_sizes_t *sizes)
{
    size_t i;

    if (sizes == NULL) {
        return;
    }

    if (sizes->names != NULL) {
        for (i = 0; i < sizes->count; i++) {
            free(sizes->names[i]);
        }
        free(sizes->names);
    }
    free(sizes->lengths);
    sizes->names = NULL;
    sizes->lengths = NULL;
    sizes->count = 0;
}

pebble_io_status_t pebble_read_chrom_sizes(const char *path, pebble_chrom_sizes_t *out)
{
    (void)path;
    if (out != NULL) {
        out->names = NULL;
        out->lengths = NULL;
        out->count = 0;
    }
    return PEBBLE_IO_ERR_UNSUPPORTED;
}

int pebble_chrom_sizes_contains(const pebble_chrom_sizes_t *sizes, const char *chrom)
{
    (void)sizes;
    (void)chrom;
    return 0;
}

pebble_output_format_t pebble_infer_output_format(const char *path)
{
    size_t len;

    if (path == NULL) {
        return PEBBLE_OUTPUT_BEDGRAPH;
    }

    len = strlen(path);
    if (len >= 3 && strcmp(path + len - 3, ".bw") == 0) {
        return PEBBLE_OUTPUT_BIGWIG;
    }
    if (len >= 8 && strcmp(path + len - 8, ".bigwig") == 0) {
        return PEBBLE_OUTPUT_BIGWIG;
    }
    return PEBBLE_OUTPUT_BEDGRAPH;
}

pebble_bigwig_writer_t *pebble_bigwig_writer_create(
    const char *path,
    const pebble_chrom_sizes_t *sizes)
{
    (void)path;
    (void)sizes;
    return NULL;
}

pebble_io_status_t pebble_bigwig_write_bedgraph(
    pebble_bigwig_writer_t *writer,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count)
{
    (void)writer;
    (void)chrom;
    (void)start_offset;
    (void)config;
    (void)values;
    (void)value_count;
    return PEBBLE_IO_ERR_UNSUPPORTED;
}

pebble_io_status_t pebble_bigwig_write_records(
    pebble_bigwig_writer_t *writer,
    const pebble_bedgraph_record_t *records,
    size_t count)
{
    (void)writer;
    (void)records;
    (void)count;
    return PEBBLE_IO_ERR_UNSUPPORTED;
}

pebble_io_status_t pebble_bigwig_writer_close(pebble_bigwig_writer_t *writer)
{
    (void)writer;
    return PEBBLE_IO_ERR_UNSUPPORTED;
}
