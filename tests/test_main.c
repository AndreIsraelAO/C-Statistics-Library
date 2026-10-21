#include <math.h>
#include <stdio.h>
#include "cstats.h"

static int g_checks = 0;
static int g_failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        g_checks++;                                                          \
        if (!(cond)) {                                                       \
            g_failures++;                                                    \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);  \
        }                                                                    \
    } while (0)

#define CHECK_NEAR(actual, expected) CHECK(fabs((actual) - (expected)) < 1e-9)
#define CHECK_NAN(actual) CHECK(isnan(actual))

/* Valores de referência gerados com o módulo `statistics` do Python. */
static const double X[] = {1, 2, 2, 3, 4, 5, 5, 5, 6};
static const double Y[] = {2, 3, 2, 5, 4, 6, 5, 7, 6};
#define N (sizeof(X) / sizeof(X[0]))

static void test_central_tendency(void) {
    CHECK_NEAR(cstats_mean(X, N), 3.6666666666666665);
    CHECK_NEAR(cstats_median(X, N), 4.0);
    CHECK_NEAR(cstats_median_low(X, N), 4.0);
    CHECK_NEAR(cstats_median_high(X, N), 4.0);
    CHECK_NEAR(cstats_mode(X, N), 5.0);
    CHECK_NEAR(cstats_harmonic_mean(X, N), 2.6865671641791047);
    CHECK_NEAR(cstats_geometric_mean(X, N), 3.208158624336986);
    CHECK_NEAR(cstats_median_grouped(X, N, 1.0), 4.0);

    const double even[] = {1, 2, 3, 4};
    CHECK_NEAR(cstats_median(even, 4), 2.5);
    CHECK_NEAR(cstats_median_low(even, 4), 2.0);
    CHECK_NEAR(cstats_median_high(even, 4), 3.0);
}

static void test_multimode(void) {
    double modes[4];
    CHECK(cstats_multimode(X, N, modes, 4) == 1);
    CHECK_NEAR(modes[0], 5.0);

    const double bi[] = {1, 1, 2, 2, 3};
    CHECK(cstats_multimode(bi, 5, modes, 4) == 2);
    CHECK_NEAR(modes[0], 1.0);
    CHECK_NEAR(modes[1], 2.0);

    /* buffer pequeno: devolve o total, grava só o que cabe */
    CHECK(cstats_multimode(bi, 5, modes, 1) == 2);
    CHECK_NEAR(modes[0], 1.0);
}

static void test_quantiles(void) {
    CHECK_NEAR(cstats_quantile(X, N, 0.25), 2.0);
    CHECK_NEAR(cstats_quantile(X, N, 0.50), 4.0);
    CHECK_NEAR(cstats_quantile(X, N, 0.75), 5.0);
    CHECK_NEAR(cstats_quantile(X, N, 0.0), 1.0);
    CHECK_NEAR(cstats_quantile(X, N, 1.0), 6.0);

    const double q[] = {0.25, 0.5, 0.75};
    double r[3];
    CHECK(cstats_quantiles(X, N, q, 3, r) == 0);
    CHECK_NEAR(r[0], 2.0);
    CHECK_NEAR(r[1], 4.0);
    CHECK_NEAR(r[2], 5.0);

    const double bad[] = {0.5, 1.5};
    CHECK(cstats_quantiles(X, N, bad, 2, r) == -1);
}

static void test_dispersion(void) {
    CHECK_NEAR(cstats_variance(X, N), 3.0);
    CHECK_NEAR(cstats_pvariance(X, N), 2.6666666666666665);
    CHECK_NEAR(cstats_stdev(X, N), 1.7320508075688772);
    CHECK_NEAR(cstats_pstdev(X, N), 1.632993161855452);
}

static void test_relations(void) {
    CHECK_NEAR(cstats_covariance(X, Y, N), 2.7916666666666665);
    CHECK_NEAR(cstats_correlation(X, Y, N), 0.8902524557833292);

    double slope, intercept, r2, resvar;
    CHECK(cstats_linear_regression(X, Y, N, &slope, &intercept) == 0);
    CHECK_NEAR(slope, 0.9305555555555555);
    CHECK_NEAR(intercept, 1.0324074074074079);

    CHECK(cstats_linear_regression_full(X, Y, N, &slope, &intercept, &r2, &resvar) == 0);
    CHECK_NEAR(r2, 0.7925494350282486);
    CHECK_NEAR(resvar, 0.7771164021164021);

    /* x constante: regressão indefinida */
    const double flat[] = {2, 2, 2};
    const double any[] = {1, 2, 3};
    CHECK(cstats_linear_regression(flat, any, 3, &slope, &intercept) == -1);
    CHECK_NAN(cstats_correlation(flat, any, 3));

    /* n == 2: reta exata, variância residual indefinida (n - 2 == 0) */
    const double x2[] = {1, 2};
    const double y2[] = {3, 5};
    CHECK(cstats_linear_regression_full(x2, y2, 2, &slope, &intercept, &r2, &resvar) == 0);
    CHECK_NEAR(slope, 2.0);
    CHECK_NAN(resvar);
}

static void test_invalid_input(void) {
    CHECK_NAN(cstats_mean(NULL, 3));
    CHECK_NAN(cstats_mean(X, 0));
    CHECK_NAN(cstats_median(NULL, 3));
    CHECK_NAN(cstats_median(X, 0));
    CHECK_NAN(cstats_variance(X, 1));
    CHECK_NAN(cstats_covariance(X, NULL, N));
    CHECK_NAN(cstats_quantile(X, N, -0.1));
    CHECK_NAN(cstats_quantile(X, N, 1.1));
    CHECK_NAN(cstats_quantile(X, N, NAN));
    CHECK_NAN(cstats_median_grouped(X, N, 0.0));
    CHECK(cstats_multimode(X, 0, NULL, 0) == -1);

    const double neg[] = {1, -2, 3};
    CHECK_NAN(cstats_harmonic_mean(neg, 3));
    CHECK_NAN(cstats_geometric_mean(neg, 3));

    CHECK_NEAR(cstats_pvariance(X, 1), 0.0);
}

static void test_nan_in_data(void) {
    /* NaN nos dados não pode gerar resultado de ordenação indefinida */
    const double d[] = {1, NAN, 3, 2};
    CHECK_NAN(cstats_median(d, 4));
    CHECK_NAN(cstats_median_low(d, 4));
    CHECK_NAN(cstats_median_high(d, 4));
    CHECK_NAN(cstats_mode(d, 4));
    CHECK_NAN(cstats_quantile(d, 4, 0.5));
    CHECK_NAN(cstats_median_grouped(d, 4, 1.0));
    CHECK_NAN(cstats_mean(d, 4));

    double modes[2];
    CHECK(cstats_multimode(d, 4, modes, 2) == -1);
}

static void test_shape(void) {
    CHECK_NEAR(cstats_skewness(X, N), -0.22113449066792545);
    CHECK_NEAR(cstats_kurtosis(X, N), -1.3541666666666667);

    const double flat[] = {4, 4, 4};
    CHECK_NAN(cstats_skewness(flat, 3));
    CHECK_NAN(cstats_kurtosis(flat, 3));
    CHECK_NAN(cstats_skewness(NULL, 3));
    CHECK_NAN(cstats_kurtosis(X, 0));
}

static void test_accumulator(void) {
    cstats_accumulator acc;
    cstats_accumulator_init(&acc);
    CHECK(cstats_accumulator_count(&acc) == 0);
    CHECK_NAN(cstats_accumulator_mean(&acc));
    CHECK_NAN(cstats_accumulator_min(&acc));

    for (size_t i = 0; i < N; i++) CHECK(cstats_accumulator_add(&acc, X[i]) == 0);
    CHECK(cstats_accumulator_count(&acc) == N);
    CHECK_NEAR(cstats_accumulator_mean(&acc), cstats_mean(X, N));
    CHECK_NEAR(cstats_accumulator_variance(&acc), cstats_variance(X, N));
    CHECK_NEAR(cstats_accumulator_pvariance(&acc), cstats_pvariance(X, N));
    CHECK_NEAR(cstats_accumulator_stdev(&acc), cstats_stdev(X, N));
    CHECK_NEAR(cstats_accumulator_min(&acc), 1.0);
    CHECK_NEAR(cstats_accumulator_max(&acc), 6.0);

    /* merge de duas metades == acumulador único */
    cstats_accumulator a, b;
    cstats_accumulator_init(&a);
    cstats_accumulator_init(&b);
    for (size_t i = 0; i < 4; i++) cstats_accumulator_add(&a, X[i]);
    for (size_t i = 4; i < N; i++) cstats_accumulator_add(&b, X[i]);
    CHECK(cstats_accumulator_merge(&a, &b) == 0);
    CHECK(cstats_accumulator_count(&a) == N);
    CHECK_NEAR(cstats_accumulator_mean(&a), cstats_mean(X, N));
    CHECK_NEAR(cstats_accumulator_variance(&a), cstats_variance(X, N));
    CHECK_NEAR(cstats_accumulator_min(&a), 1.0);
    CHECK_NEAR(cstats_accumulator_max(&a), 6.0);

    /* entradas inválidas */
    CHECK(cstats_accumulator_add(NULL, 1.0) == -1);
    CHECK(cstats_accumulator_add(&acc, NAN) == -1);
    CHECK(cstats_accumulator_count(&acc) == N);
    CHECK(cstats_accumulator_merge(NULL, &a) == -1);

    cstats_accumulator one;
    cstats_accumulator_init(&one);
    cstats_accumulator_add(&one, 7.0);
    CHECK_NAN(cstats_accumulator_variance(&one));
    CHECK_NEAR(cstats_accumulator_pvariance(&one), 0.0);
}

int main(void) {
    test_central_tendency();
    test_multimode();
    test_quantiles();
    test_dispersion();
    test_relations();
    test_invalid_input();
    test_nan_in_data();
    test_shape();
    test_accumulator();

    printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
