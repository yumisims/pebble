#ifndef PEBBLE_PEBBLE_H
#define PEBBLE_PEBBLE_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    int window_size; /* W = 1000 */
    int step_size;   /* S = 100 */
    double trim_low; /* 0.40 (trim bottom 40%) */
    double trim_high; /* 0.40 (trim top 40%) */
} pebble_config_t;

typedef enum {
    PEBBLE_OK = 0,
    PEBBLE_ERR_INVALID_ARG = -1,
    PEBBLE_ERR_TOO_SHORT = -2,
    PEBBLE_ERR_BUFFER_TOO_SMALL = -3
} pebble_status_t;

/*
 * Returns the number of smoothed values that pebble_process would write.
 * Returns 0 when coverage is shorter than the window size or config is invalid.
 */
size_t pebble_output_count(size_t coverage_len, const pebble_config_t *config);

/*
 * Processes a single chromosome's raw coverage array using a sliding trimmed mean.
 * Each window uses a coverage histogram (not sorting) to compute the trimmed mean.
 *
 * out must point to at least pebble_output_count(...) elements.
 * On success, *out_len receives the number of values written to out.
 *
 * The caller retains ownership of all buffers.
 */
pebble_status_t pebble_process(
    const int16_t *coverage,
    size_t coverage_len,
    const pebble_config_t *config,
    double *out,
    size_t out_cap,
    size_t *out_len
);

#endif /* PEBBLE_PEBBLE_H */
