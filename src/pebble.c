#include "pebble/pebble.h"

#include <stdint.h>
#include <stdlib.h>

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

static void histogram_add_range(
    const int16_t *coverage,
    size_t from,
    size_t to,
    uint32_t *counts,
    int *hist_max)
{
    size_t k;

    for (k = from; k < to; k++) {
        int idx = (int)(uint16_t)coverage[k];

        counts[idx]++;
        if (idx > *hist_max) {
            *hist_max = idx;
        }
    }
}

static void histogram_remove_range(
    const int16_t *coverage,
    size_t from,
    size_t to,
    uint32_t *counts,
    int *hist_max)
{
    size_t k;

    for (k = from; k < to; k++) {
        int idx = (int)(uint16_t)coverage[k];

        counts[idx]--;
    }

    while (*hist_max > 0 && counts[*hist_max] == 0U) {
        (*hist_max)--;
    }
}

static double trimmed_mean_from_histogram(
    const uint32_t *counts,
    int hist_max,
    int low_idx,
    int high_idx)
{
    int need_skip = low_idx;
    int need_take = high_idx - low_idx;
    int64_t sum = 0;
    int v;

    for (v = 0; v <= hist_max && need_take > 0; v++) {
        uint32_t c = counts[v];

        if (c == 0U) {
            continue;
        }

        if (need_skip > 0) {
            if ((int)c <= need_skip) {
                need_skip -= (int)c;
                continue;
            }

            {
                int avail = (int)c - need_skip;
                int take = avail < need_take ? avail : need_take;

                need_skip = 0;
                sum += (int64_t)v * take;
                need_take -= take;
            }
        } else {
            int take = (int)c < need_take ? (int)c : need_take;

            sum += (int64_t)v * take;
            need_take -= take;
        }
    }

    return (double)sum / (double)(high_idx - low_idx);
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
    uint32_t *counts = NULL;
    int hist_max = 0;
    size_t window_size;
    size_t step_size;

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

    window_size = (size_t)config->window_size;
    step_size = (size_t)config->step_size;

    counts = (uint32_t *)calloc(65536U, sizeof(uint32_t));
    if (counts == NULL) {
        return PEBBLE_ERR_INVALID_ARG;
    }

    for (size_t i = 0; i < num_steps; i++) {
        size_t start = i * step_size;

        if (i == 0U) {
            histogram_add_range(coverage, start, start + window_size, counts, &hist_max);
        } else {
            histogram_remove_range(coverage, start - step_size, start, counts, &hist_max);
            histogram_add_range(
                coverage,
                start + window_size - step_size,
                start + window_size,
                counts,
                &hist_max
            );
        }

        out[i] = trimmed_mean_from_histogram(counts, hist_max, low_idx, high_idx);
    }

    free(counts);
    *out_len = num_steps;
    return PEBBLE_OK;
}
