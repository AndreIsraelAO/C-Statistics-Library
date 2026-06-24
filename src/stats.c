/* =========================================================================
 * stats.c - Implementation of dynamic series data structure
 * =========================================================================
 * This file implements the internal structure and lifecycle functions for
 * the StatSeries opaque data type.
 */

#include "stats.h"
#include <stdlib.h>
#include <string.h>

/* =========================================================================
 * INTERNAL STRUCTURE DEFINITION
 * ========================================================================= */

/**
 * Internal definition of StatSeries - hidden from header users.
 * Contains a dynamically allocated float array with capacity tracking.
 */
struct StatSeries {
    float* data;      /* Dynamically allocated array of values */
    size_t size;      /* Current number of elements */
    size_t capacity;  /* Allocated capacity of the buffer */
};

/* =========================================================================
 * LIFECYCLE FUNCTIONS
 * ========================================================================= */

StatSeries* series_create(size_t initial_capacity) {
    /* Allocate the struct itself */
    StatSeries* s = (StatSeries*)malloc(sizeof(StatSeries));
    if (!s) {
        return NULL;
    }

    /* Handle zero initial capacity gracefully */
    if (initial_capacity == 0) {
        initial_capacity = 4; /* Default minimum capacity */
    }

    /* Allocate the internal buffer */
    s->data = (float*)malloc(initial_capacity * sizeof(float));
    if (!s->data) {
        free(s);
        return NULL;
    }

    s->size = 0;
    s->capacity = initial_capacity;

    return s;
}

void series_destroy(StatSeries* s) {
    if (!s) {
        return; /* Safely ignore NULL */
    }

    /* Free internal buffer first, then the struct */
    free(s->data);
    free(s);
}

/* =========================================================================
 * DATA MANAGEMENT FUNCTIONS
 * ========================================================================= */

bool series_add_value(StatSeries* s, float value) {
    if (!s) {
        return false;
    }

    /* Check if we need to resize */
    if (s->size >= s->capacity) {
        /* Calculate new capacity (double the current) */
        size_t new_capacity = s->capacity * 2;
        
        /* Reallocate buffer */
        float* new_data = (float*)realloc(s->data, new_capacity * sizeof(float));
        if (!new_data) {
            return false; /* Allocation failed */
        }

        s->data = new_data;
        s->capacity = new_capacity;
    }

    /* Add the value */
    s->data[s->size] = value;
    s->size++;

    return true;
}

size_t series_size(const StatSeries* s) {
    if (!s) {
        return 0;
    }
    return s->size;
}

const float* series_data(const StatSeries* s) {
    if (!s) {
        return NULL;
    }
    return s->data;
}
