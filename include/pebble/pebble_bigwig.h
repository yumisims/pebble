#ifndef PEBBLE_PEBBLE_BIGWIG_H
#define PEBBLE_PEBBLE_BIGWIG_H

#include "pebble/pebble.h"
#include "pebble/pebble_io.h"

#include <stdint.h>

typedef enum {
    PEBBLE_OUTPUT_BEDGRAPH = 0,
    PEBBLE_OUTPUT_BIGWIG = 1
} pebble_output_format_t;

typedef struct {
    char **names;
    uint32_t *lengths;
    size_t count;
} pebble_chrom_sizes_t;

void pebble_chrom_sizes_free(pebble_chrom_sizes_t *sizes);

pebble_io_status_t pebble_read_chrom_sizes(const char *path, pebble_chrom_sizes_t *out);

int pebble_chrom_sizes_contains(const pebble_chrom_sizes_t *sizes, const char *chrom);

pebble_output_format_t pebble_infer_output_format(const char *path);

typedef struct pebble_bigwig_writer pebble_bigwig_writer_t;

pebble_bigwig_writer_t *pebble_bigwig_writer_create(
    const char *path,
    const pebble_chrom_sizes_t *sizes
);

pebble_io_status_t pebble_bigwig_write_bedgraph(
    pebble_bigwig_writer_t *writer,
    const char *chrom,
    int start_offset,
    const pebble_config_t *config,
    const double *values,
    size_t value_count
);

pebble_io_status_t pebble_bigwig_write_records(
    pebble_bigwig_writer_t *writer,
    const pebble_bedgraph_record_t *records,
    size_t count
);

pebble_io_status_t pebble_bigwig_writer_close(pebble_bigwig_writer_t *writer);

#endif /* PEBBLE_PEBBLE_BIGWIG_H */
