#ifndef CSTATISTICS_STATS_H
#define CSTATISTICS_STATS_H

/**
 * @file stats.h
 * @brief Public API for dynamic series data structure and statistical computations.
 * 
 * This header provides an opaque data structure for managing dynamically allocated
 * arrays of floating-point values, along with bundled statistical result structures.
 * All functions are designed for memory safety and numerical stability.
 */

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * OPAQUE DATA STRUCTURE
 * ========================================================================= */

/**
 * @typedef StatSeries
 * @brief Opaque handle to a dynamic array of floats.
 * 
 * The internal structure is hidden from users of this API to ensure
 * encapsulation and prevent direct manipulation of internal state.
 * Use the provided lifecycle and accessor functions to interact with
 * StatSeries instances.
 */
typedef struct StatSeries StatSeries;

/* =========================================================================
 * RESULT STRUCTURE
 * ========================================================================= */

/**
 * @struct StatsResult
 * @brief Bundled statistical results from a single computation pass.
 * 
 * This structure contains commonly computed statistics grouped together
 * for convenience. All values are computed in a single pass using
 * numerically stable algorithms (Welford's method).
 * 
 * @member mean       Arithmetic mean of the data series
 * @member variance   Sample variance (using Bessel's correction, n-1)
 * @member std_dev    Sample standard deviation (sqrt of variance)
 * @member min        Minimum value in the series
 * @member max        Maximum value in the series
 * @member count      Number of elements processed
 */
typedef struct {
    float mean;       /**< Arithmetic mean of the data series */
    float variance;   /**< Sample variance (Bessel's correction applied) */
    float std_dev;    /**< Sample standard deviation (sqrt(variance)) */
    float min;        /**< Minimum value in the series */
    float max;        /**< Maximum value in the series */
    size_t count;     /**< Number of elements in the series */
} StatsResult;

/* =========================================================================
 * LIFECYCLE FUNCTIONS
 * ========================================================================= */

/**
 * @brief Creates a new StatSeries with the specified initial capacity.
 * 
 * Allocates memory for a new StatSeries instance and initializes its
 * internal buffer with the given capacity. The series starts empty
 * (size = 0) and can grow dynamically as values are added.
 * 
 * @param initial_capacity  Initial capacity of the internal buffer.
 *                          If 0, a default minimum capacity (4) is used.
 * @return Pointer to newly allocated StatSeries on success, 
 *         or NULL if memory allocation fails.
 * 
 * @note Caller is responsible for calling series_destroy() to free
 *       the allocated memory when done.
 * @note If initial_capacity is 0, it will be set to a default value (4).
 * 
 * @code
 * StatSeries* s = series_create(10);
 * if (!s) {
 *     // Handle allocation failure
 * }
 * @endcode
 */
StatSeries* series_create(size_t initial_capacity);

/**
 * @brief Destroys a StatSeries, freeing all associated memory.
 * 
 * Frees both the internal data buffer and the StatSeries structure itself.
 * After calling this function, the pointer must not be used again.
 * 
 * @param s  Pointer to StatSeries to destroy.
 *           If NULL, the function returns immediately (no-op).
 * 
 * @note This function safely handles NULL pointers.
 * @note After destruction, any pointers obtained via series_data() 
 *       become invalid.
 * 
 * @code
 * series_destroy(s);
 * s = NULL;  // Good practice to avoid dangling pointers
 * @endcode
 */
void series_destroy(StatSeries* s);

/* =========================================================================
 * DATA MANAGEMENT FUNCTIONS
 * ========================================================================= */

/**
 * @brief Appends a value to the series.
 * 
 * Adds a new floating-point value to the end of the series. If the current
 * size equals or exceeds the capacity, the internal buffer is automatically
 * resized (typically doubled) to accommodate the new element.
 * 
 * @param s      Pointer to StatSeries (must not be NULL).
 * @param value  The float value to append.
 * @return true on success, false if memory reallocation fails.
 * 
 * @note If s is NULL, returns false immediately.
 * @note On allocation failure, the series remains unchanged (no partial writes).
 * @note The internal buffer grows by doubling its capacity when full.
 * 
 * @code
 * if (!series_add_value(series, 3.14f)) {
 *     // Handle allocation failure
 * }
 * @endcode
 */
bool series_add_value(StatSeries* s, float value);

/**
 * @brief Returns the current number of elements in the series.
 * 
 * Retrieves the count of values that have been added to the series.
 * This may be less than the allocated capacity.
 * 
 * @param s  Pointer to StatSeries (must not be NULL).
 * @return Current element count, or 0 if s is NULL.
 * 
 * @note Returns 0 for NULL input (does not crash).
 * @note This is the number of valid elements, not the buffer capacity.
 * 
 * @code
 * size_t n = series_size(series);
 * printf("Series has %zu elements\n", n);
 * @endcode
 */
size_t series_size(const StatSeries* s);

/**
 * @brief Returns a read-only pointer to the internal data buffer.
 * 
 * Provides direct access to the underlying float array for efficient
 * iteration or passing to other functions. The returned pointer remains
 * valid until the series is destroyed or modified.
 * 
 * @param s  Pointer to StatSeries (must not be NULL).
 * @return Pointer to internal float array on success, or NULL if s is NULL.
 * 
 * @warning The returned pointer is READ-ONLY. Modifying the data through
 *          this pointer may lead to undefined behavior.
 * @note Returns NULL for NULL input (does not crash).
 * @note The pointer becomes invalid after series_destroy() is called.
 * 
 * @code
 * const float* data = series_data(series);
 * for (size_t i = 0; i < series_size(series); i++) {
 *     printf("%f ", data[i]);
 * }
 * @endcode
 */
const float* series_data(const StatSeries* s);

#ifdef __cplusplus
}
#endif

#endif /* CSTATISTICS_STATS_H */
