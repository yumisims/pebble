#ifndef PEBBLE_PEBBLE_IO_H
#define PEBBLE_PEBBLE_IO_H

#include "pebble/pebble.h"

#include <stdio.h>

typedef enum {
    PEBBLE_IO_OK = 0,
    PEBBLE_IO_ERR_INVALID_ARG = -1,
    PEBBLE_IO_ERR_IO = -2,
    PEBBLE_IO_ERR_PARSE = -3,
    PEBBLE_IO_ERR_OOM = -4,
    PEBBLE_IO_ERR_UNSUPPORTED = -5
} pebble_io_status_t;

typedef enum {
    PEBBLE_FORMAT_BEDGRAPH = 0,
    PEBBLE_FORMAT_BED = 1
} pebble_input_format_t;

typedef struct {
    char *chrom;
    int start_offset;
    int16_t *coverage;
    size_t length;
} pebble_coverage_t;

typedef struct {
    pebble_coverage_t *items;
    size_t count;
} pebble_coverage_batch_t;

void pebble_coverage_batch_free(pebble_coverage_batch_t *batch);

pebble_io_status_t pebble_read_bedgraph(
    const char *path,
    const char *chrom_filter,
    pebble_coverage_batch_t *out
);

pebble_io_status_t pebble_read_bed(
    const char *path,
    const char *chrom_filter,
    pebble_coverage_batch_t *out
);

pebble_io_status_t pebble_read_coverage(
    const char *path,
    pebble_input_format_t format,
    const char *chrom_filter,
    pebble_coverage_batch_t *out
);

pebble_io_status_t pebble_write_bedgraph(
    FILE *out,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count
);

pebble_io_status_t pebble_write_bedgraph_file(
    const char *path,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count
);

const char *pebble_io_status_string(pebble_io_status_t status);

#endif /* PEBBLE_PEBBLE_IO_H */
