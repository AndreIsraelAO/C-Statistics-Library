/* Micro-benchmark for cstats. Build with -DSTATS_BUILD_BENCH=ON in Release mode.
 * Usage: cstats_bench [n]   (default n = 1,000,000)
 * Prints the best of several runs per function. Timings depend heavily on the
 * machine, compiler and flags; compare only runs made under the same setup. */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cstats.h"

#define REPEATS 5

static volatile double g_sink;  /* keeps the compiler from dropping the calls */

static double now_ms(void) {
    return (double)clock() * 1000.0 / (double)CLOCKS_PER_SEC;
}

/* Deterministic pseudo-random values in [0, 1000) so runs are comparable. */
static void fill(double *a, size_t n, unsigned long long seed) {
    for (size_t i = 0; i < n; i++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        a[i] = (double)(seed >> 40) / (double)(1ULL << 24) * 1000.0;
    }
}

#define BENCH(label, count, expr)                                        \
    do {                                                                 \
        double best = 1e300;                                             \
        for (int r = 0; r < REPEATS; r++) {                              \
            double t0 = now_ms();                                        \
            expr;                                                        \
            double dt = now_ms() - t0;                                   \
            if (dt < best) best = dt;                                    \
        }                                                                \
        printf("%-28s n=%-9zu %9.3f ms\n", label, (size_t)(count), best); \
    } while (0)

int main(int argc, char **argv) {
    size_t n = argc > 1 ? (size_t)strtoull(argv[1], NULL, 10) : 1000000;
    if (n < 2) n = 2;

    double *x = (double*)malloc(n * sizeof(double));
    double *y = (double*)malloc(n * sizeof(double));
    if (!x || !y) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    fill(x, n, 1);
    fill(y, n, 2);

    size_t m = n / 10 < 2 ? 2 : n / 10;  /* Spearman sorts twice: use a smaller n */

    BENCH("mean", n, g_sink = cstats_mean(x, n));
    BENCH("variance", n, g_sink = cstats_variance(x, n));
    BENCH("correlation (Pearson)", n, g_sink = cstats_correlation(x, y, n));
    BENCH("median (sorts a copy)", n, g_sink = cstats_median(x, n));
    BENCH("quantile 0.9 (sorts a copy)", n, g_sink = cstats_quantile(x, n, 0.9));
    BENCH("mad", n, g_sink = cstats_mad(x, n));
    BENCH("spearman", m, g_sink = cstats_spearman(x, y, m));
    BENCH("accumulator add + stdev", n, {
        cstats_accumulator acc;
        cstats_accumulator_init(&acc);
        for (size_t i = 0; i < n; i++) cstats_accumulator_add(&acc, x[i]);
        g_sink = cstats_accumulator_stdev(&acc);
    });

    free(x);
    free(y);
    return 0;
}
