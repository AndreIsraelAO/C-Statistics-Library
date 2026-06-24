/**
 * example_basic.c - Basic demonstration of the C Statistics Library
 * ------------------------------------------------------------------
 * This example demonstrates how to use the StatSeries data structure
 * and the statistical functions provided by the library.
 */

#include <stdio.h>
#include "stats.h"
#include "cstats.h"

int main(void) {
    printf("=== C Statistics Library Demo ===\n\n");

    /* Create a new series with initial capacity of 10 */
    StatSeries* series = series_create(10);
    if (!series) {
        fprintf(stderr, "Failed to create series\n");
        return 1;
    }

    /* Add some sample data */
    printf("Adding values: 1.0, 2.0, 3.0, 4.0, 5.0\n");
    series_add_value(series, 1.0f);
    series_add_value(series, 2.0f);
    series_add_value(series, 3.0f);
    series_add_value(series, 4.0f);
    series_add_value(series, 5.0f);

    /* Get series info */
    printf("Series size: %zu\n", series_size(series));

    /* Calculate basic statistics using the raw data */
    const float* data = series_data(series);
    size_t count = series_size(series);

    /* Convert to double for cstats functions */
    double double_data[5];
    for (size_t i = 0; i < count; i++) {
        double_data[i] = (double)data[i];
    }

    double mean = cstats_mean(double_data, count);
    double median = cstats_median(double_data, count);
    double variance = cstats_variance(double_data, count);
    double std_dev = cstats_stdev(double_data, count);

    printf("\nStatistics:\n");
    printf("  Mean:     %.2f\n", mean);
    printf("  Median:   %.2f\n", median);
    printf("  Variance: %.2f\n", variance);
    printf("  Std Dev:  %.2f\n", std_dev);

    /* Clean up */
    series_destroy(series);

    printf("\nDemo completed successfully!\n");
    return 0;
}
