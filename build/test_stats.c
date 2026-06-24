#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "stats.h"

int main(void) {
    printf("Testing StatSeries implementation...\n");
    
    /* Test 1: Create with initial capacity */
    StatSeries* s = series_create(4);
    assert(s != NULL);
    assert(series_size(s) == 0);
    printf("PASS: series_create(4)\n");
    
    /* Test 2: Add values */
    assert(series_add_value(s, 1.0f) == true);
    assert(series_add_value(s, 2.0f) == true);
    assert(series_add_value(s, 3.0f) == true);
    assert(series_add_value(s, 4.0f) == true);
    assert(series_size(s) == 4);
    printf("PASS: series_add_value (4 values)\n");
    
    /* Test 3: Automatic resize */
    assert(series_add_value(s, 5.0f) == true);
    assert(series_add_value(s, 6.0f) == true);
    assert(series_size(s) == 6);
    printf("PASS: automatic resize\n");
    
    /* Test 4: Data access */
    const float* data = series_data(s);
    assert(data != NULL);
    assert(data[0] == 1.0f);
    assert(data[5] == 6.0f);
    printf("PASS: series_data access\n");
    
    /* Test 5: NULL handling */
    assert(series_size(NULL) == 0);
    assert(series_data(NULL) == NULL);
    assert(series_add_value(NULL, 1.0f) == false);
    series_destroy(NULL); /* Should not crash */
    printf("PASS: NULL handling\n");
    
    /* Test 6: Destroy */
    series_destroy(s);
    printf("PASS: series_destroy\n");
    
    /* Test 7: Create with zero capacity */
    s = series_create(0);
    assert(s != NULL);
    assert(series_add_value(s, 1.0f) == true);
    series_destroy(s);
    printf("PASS: series_create(0)\n");
    
    printf("\nAll tests passed!\n");
    return 0;
}
