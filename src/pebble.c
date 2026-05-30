#include "pebble/pebble.h"

#include <stdlib.h>

typedef struct {
    const int16_t *values;
} sort_ctx_t;

static int compare_int16(const void *left, const void *right)
{
    int16_t a = *(const int16_t *)left;
    int16_t b = *(const int16_t *)right;

    if (a < b) {
        return -1;
    }
    if (a > b) {
        return 1;
    }
    return 0;
}

static int config_is_valid(const pebble_config_t *config)
{
    if (config == NULL) {
        return 0;
    }
    if (config->window_size <= 0 || config->step_size <= 0) {
        return 0;
    }
    if (config->trim_low < 0.0 || config->trim_high < 0.0) {
        return 0;
    }
    if (config->trim_low + config->trim_high >= 1.0) {
        return 0;
    }
    return 1;
}

size_t pebble_output_count(size_t coverage_len, const pebble_config_t *config)
{
    if (!config_is_valid(config)) {
        return 0;
    }
    if (coverage_len < (size_t)config->window_size) {
        return 0;
    }

    return ((coverage_len - (size_t)config->window_size) / (size_t)config->step_size) + 1U;
}

pebble_status_t pebble_process(
    const int16_t *coverage,
    size_t coverage_len,
    const pebble_config_t *config,
    double *out,
    size_t out_cap,
    size_t *out_len)
{
    size_t num_steps;
    int low_idx;
    int high_idx;
    double trimmed_len;
    int16_t *window_buf = NULL;

    if (coverage == NULL || config == NULL || out == NULL || out_len == NULL) {
        return PEBBLE_ERR_INVALID_ARG;
    }

    if (!config_is_valid(config)) {
        return PEBBLE_ERR_INVALID_ARG;
    }

    if (coverage_len < (size_t)config->window_size) {
        return PEBBLE_ERR_TOO_SHORT;
    }

    num_steps = pebble_output_count(coverage_len, config);
    if (num_steps == 0U) {
        return PEBBLE_ERR_TOO_SHORT;
    }
    if (out_cap < num_steps) {
        return PEBBLE_ERR_BUFFER_TOO_SMALL;
    }

    low_idx = (int)((double)config->window_size * config->trim_low);
    high_idx = (int)((double)config->window_size * (1.0 - config->trim_high));
    trimmed_len = (double)(high_idx - low_idx);
    if (trimmed_len <= 0.0) {
        return PEBBLE_ERR_INVALID_ARG;
    }

    window_buf = (int16_t *)malloc((size_t)config->window_size * sizeof(int16_t));
    if (window_buf == NULL) {
        return PEBBLE_ERR_INVALID_ARG;
    }

    for (size_t i = 0; i < num_steps; i++) {
        size_t start = i * (size_t)config->step_size;
        int64_t sum = 0;

        for (int k = 0; k < config->window_size; k++) {
            window_buf[k] = coverage[start + (size_t)k];
        }

        qsort(window_buf, (size_t)config->window_size, sizeof(int16_t), compare_int16);

        for (int k = low_idx; k < high_idx; k++) {
            sum += (int64_t)window_buf[k];
        }

        out[i] = (double)sum / trimmed_len;
    }

    free(window_buf);
    *out_len = num_steps;
    return PEBBLE_OK;
}
