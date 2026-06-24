#ifndef CSTATISTICS_STATISTICS_H
#define CSTATISTICS_STATISTICS_H

/**
 * @file cstats.h
 * @brief Public API for statistical computations in C.
 * 
 * This header defines the stable public interface of the C Statistics Library.
 * Functions are designed to be simple, reliable, deterministic, and compatible
 * with both C and C++. All functions handle edge cases gracefully, returning
 * NaN or error codes for invalid inputs.
 */

#include <stddef.h>  /**< for size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * BASIC METRICS
 * ========================================================================= */

/**
 * @brief Calculates the arithmetic mean of a data series.
 * 
 * Computes the sum of all elements divided by the count.
 * Uses a single pass through the data for efficiency.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Arithmetic mean on success, or NaN if size == 0 or data is NULL.
 * 
 * @note Time complexity: O(n)
 * @note Returns NaN for empty arrays (size == 0).
 */
double cstats_mean(const double *data, size_t size);

/**
 * @brief Calculates the mode (most frequent value) of a data series.
 * 
 * Finds the value that appears most frequently in the dataset.
 * If multiple values share the same maximum frequency, returns
 * the smallest one (in ascending order).
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Mode value on success, or NaN if size == 0 or data is NULL.
 * 
 * @note Time complexity: O(n log n) due to sorting.
 * @note Creates a temporary copy of the data to preserve original order.
 */
double cstats_mode(const double *data, size_t size);

/**
 * @brief Calculates all modes (multimodal values) of a data series.
 * 
 * Identifies all values that appear with the maximum frequency.
 * Useful for detecting multimodal distributions.
 * 
 * @param data       Pointer to array of double values (must not be NULL).
 * @param size       Number of elements in the array.
 * @param modes      Output array to receive the mode values (must not be NULL).
 * @param max_modes  Maximum number of modes that can be stored in modes array.
 * @return Number of modes found on success, or -1 on error.
 * 
 * @note Time complexity: O(n log n) due to sorting.
 * @note Caller must allocate the modes array with sufficient capacity.
 */
int cstats_multimode(const double *data, size_t size, double *modes, size_t max_modes);

/**
 * @brief Calculates the median of a data series.
 * 
 * Finds the middle value when data is sorted. For even-sized datasets,
 * returns the average of the two central values.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Median value on success, or NaN if size == 0 or data is NULL.
 * 
 * @note Time complexity: O(n log n) due to sorting.
 * @note Creates a temporary copy to preserve original data order.
 */
double cstats_median(const double *data, size_t size);


/* =========================================================================
 * VARIOUS MEANS
 * ========================================================================= */

/**
 * @brief Calculates the harmonic mean of a data series.
 * 
 * Computes n / (sum 1/x_i), useful for rates and ratios.
 * Requires all values to be strictly positive.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Harmonic mean on success, or NaN if any value <= 0.
 * 
 * @note Time complexity: O(n)
 * @note All values must be > 0; returns NaN otherwise.
 */
double cstats_harmonic_mean(const double *data, size_t size);

/**
 * @brief Calculates the geometric mean of a data series.
 * 
 * Computes (product x_i)^(1/n), using log-sum for numerical stability.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Geometric mean on success, or NaN if any value <= 0.
 * 
 * @note Time complexity: O(n)
 * @note Uses log-sum method to avoid overflow.
 */
double cstats_geometric_mean(const double *data, size_t size);


/* =========================================================================
 * SPECIALIZED MEDIANS
 * ========================================================================= */

/**
 * @brief Calculates the low median of a data series.
 * 
 * For even-sized datasets, returns the lower of the two central values.
 * Example: [1, 2, 3, 4] -> median_low = 2
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Low median on success, or NaN if size == 0 or data is NULL.
 */
double cstats_median_low(const double *data, size_t size);

/**
 * @brief Calculates the high median of a data series.
 * 
 * For even-sized datasets, returns the higher of the two central values.
 * Example: [1, 2, 3, 4] -> median_high = 3
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return High median on success, or NaN if size == 0 or data is NULL.
 */
double cstats_median_high(const double *data, size_t size);

/**
 * @brief Calculates the grouped median of a data series.
 * 
 * Computes the median for continuous data grouped into intervals.
 * Formula: median = L + h * ((n/2 - cf) / f)
 * 
 * @param data     Pointer to array of double values (must not be NULL).
 * @param size     Number of elements in the array.
 * @param interval Width of each grouping interval (must be > 0).
 * @return Grouped median on success, or NaN if invalid input.
 * 
 * @note Mimics Python's statistics.median_grouped behavior.
 */
double cstats_median_grouped(const double *data, size_t size, double interval);

/* =========================================================================
 * QUANTILES
 * ========================================================================= */

/**
 * @brief Calculates an arbitrary quantile of a data series.
 * 
 * Returns the value below which a given fraction of observations fall.
 * Uses linear interpolation between adjacent values when needed.
 * 
 * @param data     Pointer to array of double values (must not be NULL).
 * @param size     Number of elements in the array.
 * @param quantile Quantile value in range [0.0, 1.0].
 * @return Quantile value on success, or NaN if invalid input.
 * 
 * @note Common quantiles: Q1=0.25, Q2=0.5 (median), Q3=0.75.
 */
double cstats_quantile(const double *data, size_t size, double quantile);

/**
 * @brief Calculates multiple quantiles in a single call.
 * 
 * Efficiently computes several quantiles at once, avoiding repeated sorting.
 * 
 * @param data    Pointer to array of double values (must not be NULL).
 * @param size    Number of elements in the array.
 * @param q       Array of quantile values in range [0.0, 1.0].
 * @param q_count Number of quantiles to compute.
 * @param result  Output array to receive quantile values.
 * @return 0 on success, -1 on error.
 */
int cstats_quantiles(const double *data, size_t size,
                     const double *q, size_t q_count, double *result);


/* =========================================================================
 * VARIANCE AND STANDARD DEVIATION
 * ========================================================================= */

/**
 * @brief Calculates the sample variance of a data series.
 * 
 * Computes sum((x - mean)^2) / (n - 1), using Bessel's correction.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Sample variance on success, or NaN if size < 2.
 * 
 * @note Requires at least 2 elements.
 * @note Uses Bessel's correction (divides by n-1).
 */
double cstats_variance(const double *data, size_t size);

/**
 * @brief Calculates the population variance of a data series.
 * 
 * Computes sum((x - mean)^2) / n.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Population variance on success, or NaN if size < 1.
 */
double cstats_pvariance(const double *data, size_t size);

/**
 * @brief Calculates the sample standard deviation.
 * 
 * Square root of the sample variance.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Sample standard deviation on success, or NaN if size < 2.
 */
double cstats_stdev(const double *data, size_t size);

/**
 * @brief Calculates the population standard deviation.
 * 
 * Square root of the population variance.
 * 
 * @param data  Pointer to array of double values (must not be NULL).
 * @param size  Number of elements in the array.
 * @return Population standard deviation on success, or NaN if size < 1.
 */
double cstats_pstdev(const double *data, size_t size);


/* =========================================================================
 * COVARIANCE AND CORRELATION
 * ========================================================================= */

/**
 * @brief Calculates the sample covariance between two data series.
 * 
 * Computes sum((xi - mean_x)(yi - mean_y)) / (n - 1).
 * 
 * @param x     Pointer to first array (must not be NULL).
 * @param y     Pointer to second array (must not be NULL).
 * @param size  Number of elements in each array.
 * @return Sample covariance on success, or NaN if size < 2.
 */
double cstats_covariance(const double *x, const double *y, size_t size);

/**
 * @brief Calculates Pearson's correlation coefficient.
 * 
 * Computes cov(x,y) / (stdev(x) * stdev(y)).
 * Normalized measure in range [-1, 1].
 * 
 * @param x     Pointer to first array (must not be NULL).
 * @param y     Pointer to second array (must not be NULL).
 * @param size  Number of elements in each array.
 * @return Correlation coefficient in [-1, 1], or NaN if invalid.
 * 
 * @note Returns NaN if either variable has zero variance.
 */
double cstats_correlation(const double *x, const double *y, size_t size);


/* =========================================================================
 * LINEAR REGRESSION
 * ========================================================================= */

/**
 * @brief Performs simple linear regression (y = slope * x + intercept).
 * 
 * Fits a straight line using least squares method.
 * 
 * @param x          Pointer to independent variable array.
 * @param y          Pointer to dependent variable array.
 * @param size       Number of elements in each array.
 * @param slope      Output: slope of regression line.
 * @param intercept  Output: y-intercept.
 * @return 0 on success, -1 on error (size < 2 or NULL pointers).
 */
int cstats_linear_regression(const double *x, const double *y, size_t size,
                             double *slope, double *intercept);

/**
 * @brief Performs linear regression with goodness-of-fit metrics.
 * 
 * Extended version computing R-squared and residual variance.
 * 
 * @param x                 Pointer to independent variable array.
 * @param y                 Pointer to dependent variable array.
 * @param size              Number of elements in each array.
 * @param slope             Output: slope (can be NULL).
 * @param intercept         Output: y-intercept (can be NULL).
 * @param r_squared         Output: R-squared coefficient (can be NULL).
 * @param residual_variance Output: variance of residuals (can be NULL).
 * @return 0 on success, -1 on error.
 * 
 * @note R-squared ranges from 0 to 1; higher values indicate better fit.
 */
int cstats_linear_regression_full(const double *x, const double *y, size_t size,
                                  double *slope, double *intercept,
                                  double *r_squared, double *residual_variance);


#ifdef __cplusplus
}
#endif

#endif /* CSTATISTICS_STATISTICS_H */
