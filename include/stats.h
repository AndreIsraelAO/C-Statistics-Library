#ifndef CSTATISTICS_STATS_H
#define CSTATISTICS_STATS_H

/**
 * stats.h - Public API for dynamic series data structure
 * -------------------------------------------------------
 * This header provides an opaque data structure for managing
 * dynamically allocated arrays of floating-point values.
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
 * StatSeries - Opaque handle to a dynamic array of floats.
 * The internal structure is hidden from users of this API.
 */
typedef struct StatSeries StatSeries;

/* =========================================================================
 * LIFECYCLE FUNCTIONS
 * ========================================================================= */

/**
 * Creates a new StatSeries with the specified initial capacity.
 * 
 * @param initial_capacity  Initial capacity of the internal buffer
 * @return Pointer to new StatSeries, or NULL on allocation failure
 */
StatSeries* series_create(size_t initial_capacity);

/**
 * Destroys a StatSeries, freeing all associated memory.
 * 
 * @param s  Pointer to StatSeries to destroy (NULL is safely ignored)
 */
void series_destroy(StatSeries* s);

/* =========================================================================
 * DATA MANAGEMENT FUNCTIONS
 * ========================================================================= */

/**
 * Appends a value to the series.
 * Automatically resizes the internal buffer if capacity is reached.
 * 
 * @param s      Pointer to StatSeries
 * @param value  Value to append
 * @return true on success, false on allocation failure
 */
bool series_add_value(StatSeries* s, float value);

/**
 * Returns the current number of elements in the series.
 * 
 * @param s  Pointer to StatSeries (must not be NULL)
 * @return Current element count, or 0 if s is NULL
 */
size_t series_size(const StatSeries* s);

/**
 * Returns a read-only pointer to the internal data buffer.
 * 
 * @param s  Pointer to StatSeries (must not be NULL)
 * @return Pointer to internal float array, or NULL if s is NULL
 */
const float* series_data(const StatSeries* s);

#ifdef __cplusplus
}
#endif

#endif /* CSTATISTICS_STATS_H */
