#include <stdio.h>
#include <math.h>
#include "cstats.h"

int main(void) {
    double data[] = {1, 2, 2, 3, 4, 5, 5, 5, 6};
    size_t size = sizeof(data) / sizeof(data[0]);

    printf("=== Teste de estatísticas básicas ===\n");
    printf("Mean: %f\n", cstats_mean(data, size));
    printf("Median: %f\n", cstats_median(data, size));
    printf("Median low: %f\n", cstats_median_low(data, size));
    printf("Median high: %f\n", cstats_median_high(data, size));

    double modes[10];
    int n_modes = cstats_multimode(data, size, modes, 10);
    printf("Multimode (%d): ", n_modes);
    for (int i = 0; i < n_modes; i++) printf("%f ", modes[i]);
    printf("\n");
    printf("Mode: %f\n", cstats_mode(data, size));

    printf("\n=== Teste de médias especiais ===\n");
    printf("Harmonic mean: %f\n", cstats_harmonic_mean(data, size));
    printf("Geometric mean: %f\n", cstats_geometric_mean(data, size));

    printf("\n=== Teste de medianas agrupadas ===\n");
    printf("Median grouped (interval=1.0): %f\n", cstats_median_grouped(data, size, 1.0));

    printf("\n=== Teste de quantis ===\n");
    printf("Quantile 0.25: %f\n", cstats_quantile(data, size, 0.25));
    printf("Quantile 0.5: %f\n", cstats_quantile(data, size, 0.5));
    printf("Quantile 0.75: %f\n", cstats_quantile(data, size, 0.75));

    double q_values[] = {0.25, 0.5, 0.75};
    double q_results[3];
    cstats_quantiles(data, size, q_values, 3, q_results);
    printf("Quantiles multiple: ");
    for (int i = 0; i < 3; i++) printf("%f ", q_results[i]);
    printf("\n");

    printf("\n=== Teste de variância e desvio padrão ===\n");
    printf("Variance (sample): %f\n", cstats_variance(data, size));
    printf("Variance (pop): %f\n", cstats_pvariance(data, size));
    printf("Stdev (sample): %f\n", cstats_stdev(data, size));
    printf("Stdev (pop): %f\n", cstats_pstdev(data, size));

    printf("\n=== Teste de covariância e correlação ===\n");
    double y_data[] = {2, 3, 2, 5, 4, 6, 5, 7, 6};
    printf("Covariance: %f\n", cstats_covariance(data, y_data, size));
    printf("Correlation: %f\n", cstats_correlation(data, y_data, size));

    printf("\n=== Teste de regressão linear ===\n");
    double slope, intercept;
    if (cstats_linear_regression(data, y_data, size, &slope, &intercept) == 0) {
        printf("Linear regression: y = %fx + %f\n", slope, intercept);
    }

    double r2, residual_var;
    if (cstats_linear_regression_full(data, y_data, size,
                                      &slope, &intercept,
                                      &r2, &residual_var) == 0) {
        printf("Full regression:\n");
        printf("  y = %fx + %f\n", slope, intercept);
        printf("  R^2 = %f\n", r2);
        printf("  Residual variance = %f\n", residual_var);
    }

    return 0;
}
