#include "pebble/pebble.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int double_eq(double left, double right)
{
    return fabs(left - right) < 1e-9;
}

static void test_uniform_coverage(void)
{
    pebble_config_t config = {
        .window_size = 10,
        .step_size = 5,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };

    int16_t coverage[20];
    double out[4];
    size_t out_len = 0;

    for (size_t i = 0; i < 20; i++) {
        coverage[i] = 30;
    }

    assert(pebble_output_count(20, &config) == 3U);
    assert(pebble_process(coverage, 20, &config, out, 3, &out_len) == PEBBLE_OK);
    assert(out_len == 3U);
    assert(double_eq(out[0], 30.0));
    assert(double_eq(out[1], 30.0));
    assert(double_eq(out[2], 30.0));
}

static void test_too_short_coverage(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };

    int16_t coverage[500];
    double out[1];
    size_t out_len = 0;

    assert(pebble_output_count(500, &config) == 0U);
    assert(pebble_process(coverage, 500, &config, out, 1, &out_len) == PEBBLE_ERR_TOO_SHORT);
}

static void test_mock_scaffold_matches_go_defaults(void)
{
    pebble_config_t config = {
        .window_size = 1000,
        .step_size = 100,
        .trim_low = 0.40,
        .trim_high = 0.40,
    };

    int16_t coverage[1500];
    double out[6];
    size_t out_len = 0;

    for (size_t i = 0; i < 1500; i++) {
        coverage[i] = 30;
    }
    for (size_t i = 200; i < 250; i++) {
        coverage[i] = 0;
    }
    for (size_t i = 700; i < 850; i++) {
        coverage[i] = 999;
    }

    assert(pebble_output_count(1500, &config) == 6U);
    assert(pebble_process(coverage, 1500, &config, out, 6, &out_len) == PEBBLE_OK);
    assert(out_len == 6U);

    /* Dropout and spike are trimmed out of the central 20% mean. */
    for (size_t i = 0; i < out_len; i++) {
        assert(double_eq(out[i], 30.0));
    }
}

int main(void)
{
    test_uniform_coverage();
    test_too_short_coverage();
    test_mock_scaffold_matches_go_defaults();

    puts("all tests passed");
    return EXIT_SUCCESS;
}
