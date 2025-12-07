#include "cstats.h"
#include<stddef.h>
#include<math.h>
#include <float.h>   
#include <stdlib.h> 

/* =========================================================================
 * MÉTRICAS BÁSICAS
 * ========================================================================= */


double cstats_mean(const double *data, size_t size) {
    if (!data || size == 0) {
        // Retorna NaN em caso de entrada inválida
        return NAN;
    }

    double sum = 0.0;

    // Somatório
    for (size_t i = 0; i < size; i++) {
        sum += data[i];
    }

    // Média
    return sum / (double)size;
}

// Função de comparação para qsort (ordenação crescente)
static int compare_doubles(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

//mode
double cstats_mode(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;
    }

    // Criar uma cópia do array para ordenar
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) {
        return NAN; // falha de alocação
    }

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    qsort(sorted, size, sizeof(double), compare_doubles);

    // Variáveis para contagem
    double mode = sorted[0];
    size_t max_count = 1;
    size_t current_count = 1;

    for (size_t i = 1; i < size; i++) {
        if (sorted[i] == sorted[i - 1]) {
            current_count++;
        } else {
            current_count = 1;
        }

        if (current_count > max_count) {
            max_count = current_count;
            mode = sorted[i];
        }
    }

    free(sorted);
    return mode;
}

//multimode

int cstats_multimode(const double *data, size_t size, double *modes, size_t max_modes) {
    if (!data || size == 0 || !modes || max_modes == 0) {
        return -1; // Entrada inválida
    }

    // Criar uma cópia ordenada
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) return -1; // falha de alocação

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    qsort(sorted, size, sizeof(double), compare_doubles);

    // Contagem das frequências
    size_t max_count = 1;
    size_t current_count = 1;

    // Primeiro pass: determinar a frequência máxima
    for (size_t i = 1; i < size; i++) {
        if (sorted[i] == sorted[i - 1]) {
            current_count++;
        } else {
            if (current_count > max_count) {
                max_count = current_count;
            }
            current_count = 1;
        }
    }
    if (current_count > max_count) max_count = current_count;

    // Segundo pass: coletar todas as modas
    size_t mode_index = 0;
    current_count = 1;
    for (size_t i = 1; i <= size; i++) {
        if (i < size && sorted[i] == sorted[i - 1]) {
            current_count++;
        } else {
            if (current_count == max_count) {
                if (mode_index < max_modes) {
                    modes[mode_index++] = sorted[i - 1];
                } else {
                    break; // Não cabe mais no array modes
                }
            }
            current_count = 1;
        }
    }

    free(sorted);
    return (int)mode_index; // número de modas encontradas
}

//median

double cstats_median(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;  // entrada inválida
    }

    // Criar cópia para não modificar os dados originais
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) {
        return NAN;  // falha de alocação
    }

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    // Ordenar os dados
    qsort(sorted, size, sizeof(double), compare_doubles);

    double median;
    if (size % 2 == 1) {
        // Número ímpar de elementos
        median = sorted[size / 2];
    } else {
        // Número par de elementos → média dos dois centrais
        median = (sorted[(size / 2) - 1] + sorted[size / 2]) / 2.0;
    }

    free(sorted);
    return median;
}

/* =========================================================================
 * MÉDIAS VARIADAS
 * ========================================================================= */


//harmonic_mean

double cstats_harmonic_mean(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;  // entrada inválida
    }

    double sum_reciprocal = 0.0;

    for (size_t i = 0; i < size; i++) {
        if (data[i] <= 0.0) {
            return NAN;  // todos os elementos devem ser positivos
        }
        sum_reciprocal += 1.0 / data[i];
    }

    return (double)size / sum_reciprocal;
}

double cstats_geometric_mean(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;  // entrada inválida
    }

    double log_sum = 0.0;

    for (size_t i = 0; i < size; i++) {
        if (data[i] <= 0.0) {
            return NAN;  // todos os elementos devem ser positivos
        }
        log_sum += log(data[i]);
    }

    // Retorna a média geométrica usando exp para evitar overflow
    return exp(log_sum / (double)size);
}

/* =========================================================================
 * MEDIANAS ESPECIALIZADAS
 * ========================================================================= */


//Mediana baixa

double cstats_median_low(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN; // entrada inválida
    }

    // Criar cópia para não modificar os dados originais
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) {
        return NAN; // falha de alocação
    }

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    // Ordenar os dados
    qsort(sorted, size, sizeof(double), compare_doubles);

    double median;
    if (size % 2 == 1) {
        // Número ímpar de elementos → elemento central
        median = sorted[size / 2];
    } else {
        // Número par → retorna o elemento inferior do par central
        median = sorted[(size / 2) - 1];
    }

    free(sorted);
    return median;
}

//mediana alta

double cstats_median_high(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN; // entrada inválida
    }

    // Criar cópia para não modificar os dados originais
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) {
        return NAN; // falha de alocação
    }

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    // Ordenar os dados
    qsort(sorted, size, sizeof(double), compare_doubles);

    double median;
    if (size % 2 == 1) {
        // Número ímpar de elementos → elemento central
        median = sorted[size / 2];
    } else {
        // Número par → retorna o elemento superior do par central
        median = sorted[size / 2];
    }

    free(sorted);
    return median;
}

//mediana agrupada

double cstats_median_grouped(const double *data, size_t size, double interval) {
    if (!data || size == 0 || interval <= 0.0) {
        return NAN;  // entrada inválida
    }

    // Criar cópia para ordenar
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) return NAN;  // falha de alocação

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    qsort(sorted, size, sizeof(double), compare_doubles);

    // Determinar os limites dos grupos
    double min_val = sorted[0];
    double max_val = sorted[size - 1];
    size_t n_groups = (size_t)ceil((max_val - min_val) / interval) + 1;

    // Frequências dos grupos
    size_t *freq = (size_t*)calloc(n_groups, sizeof(size_t));
    if (!freq) {
        free(sorted);
        return NAN;
    }

    for (size_t i = 0; i < size; i++) {
        size_t index = (size_t)floor((sorted[i] - min_val) / interval);
        if (index >= n_groups) index = n_groups - 1; // segurança
        freq[index]++;
    }

    // Determinar o grupo da mediana
    size_t cumulative = 0;
    size_t median_group = 0;
    size_t n_half = size / 2;
    for (size_t i = 0; i < n_groups; i++) {
        cumulative += freq[i];
        if (cumulative > n_half) {
            median_group = i;
            break;
        }
    }

    // Frequência acumulada antes do grupo da mediana
    size_t cf = 0;
    for (size_t i = 0; i < median_group; i++) {
        cf += freq[i];
    }

    // L = limite inferior do grupo da mediana
    double L = min_val + median_group * interval;
    double f = (double)freq[median_group];
    double h = interval;

    double median = L + h * ((size / 2.0 - cf) / f);

    free(freq);
    free(sorted);

    return median;
}

/* =========================================================================
 * QUANTIS
 * ========================================================================= */

//quantil 

double cstats_quantile(const double *data, size_t size, double quantile) {
    if (!data || size == 0 || quantile < 0.0 || quantile > 1.0) {
        return NAN;  // entrada inválida
    }

    // Criar cópia para ordenar
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) return NAN;  // falha de alocação

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    qsort(sorted, size, sizeof(double), compare_doubles);

    double pos = quantile * (size - 1);
    size_t index_lower = (size_t)floor(pos);
    size_t index_upper = (size_t)ceil(pos);
    double fraction = pos - index_lower;

    double result;
    if (index_lower == index_upper) {
        result = sorted[index_lower];
    } else {
        result = sorted[index_lower] + fraction * (sorted[index_upper] - sorted[index_lower]);
    }

    free(sorted);
    return result;
}
 
//quantis

int cstats_quantiles(const double *data, size_t size,
                     const double *q, size_t q_count,
                     double *result) {
    if (!data || size == 0 || !q || q_count == 0 || !result) {
        return -1;  // entrada inválida
    }

    // Criar cópia ordenada dos dados
    double *sorted = (double*)malloc(size * sizeof(double));
    if (!sorted) return -1;  // falha de alocação

    for (size_t i = 0; i < size; i++) {
        sorted[i] = data[i];
    }

    qsort(sorted, size, sizeof(double), compare_doubles);

    // Calcular cada quantil
    for (size_t i = 0; i < q_count; i++) {
        double qi = q[i];
        if (qi < 0.0 || qi > 1.0) {
            free(sorted);
            return -1;  // quantil fora do intervalo
        }

        double pos = qi * (size - 1);
        size_t index_lower = (size_t)floor(pos);
        size_t index_upper = (size_t)ceil(pos);
        double fraction = pos - index_lower;

        if (index_lower == index_upper) {
            result[i] = sorted[index_lower];
        } else {
            result[i] = sorted[index_lower] + fraction * (sorted[index_upper] - sorted[index_lower]);
        }
    }

    free(sorted);
    return 0;  // sucesso
}

/* =========================================================================
 * VARIÂNCIA E DESVIO PADRÃO
 * ========================================================================= */

//Variância amostral 

 double cstats_variance(const double *data, size_t size) {
    if (!data || size < 2) {
        return NAN;  // entrada inválida
    }

    // Primeiro, calcular a média
    double sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum += data[i];
    }
    double mean = sum / (double)size;

    // Calcular soma dos quadrados das diferenças em relação à média
    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        sq_diff_sum += diff * diff;
    }

    return sq_diff_sum / (double)(size - 1);
}

//variancia populacional

double cstats_pvariance(const double *data, size_t size) {
    if (!data || size < 1) {
        return NAN;  // entrada inválida
    }

    // Calcular a média
    double sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum += data[i];
    }
    double mean = sum / (double)size;

    // Calcular soma dos quadrados das diferenças em relação à média
    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        sq_diff_sum += diff * diff;
    }

    return sq_diff_sum / (double)size;
}

//desvio padrao amostral

double cstats_stdev(const double *data, size_t size) {
    if (!data || size < 2) {
        return NAN;  // entrada inválida
    }

    // Calcular a variância amostral
    double sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum += data[i];
    }
    double mean = sum / (double)size;

    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        sq_diff_sum += diff * diff;
    }

    return sqrt(sq_diff_sum / (double)(size - 1));
}

//desvio padrao populacional

double cstats_pstdev(const double *data, size_t size) {
    if (!data || size < 1) {
        return NAN;  // entrada inválida
    }

    // Calcular a média
    double sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum += data[i];
    }
    double mean = sum / (double)size;

    // Calcular soma dos quadrados das diferenças em relação à média
    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        sq_diff_sum += diff * diff;
    }

    return sqrt(sq_diff_sum / (double)size);
}

/* =========================================================================
 * COVARIÂNCIA E CORRELAÇÃO
 * ========================================================================= */

 //Covariância amostral entre dois vetores

 double cstats_covariance(const double *x, const double *y, size_t size) {
    if (!x || !y || size < 2) {
        return NAN;  // entrada inválida
    }

    // Calcular médias
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / (double)size;
    double mean_y = sum_y / (double)size;

    // Calcular soma dos produtos das diferenças
    double cov_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        cov_sum += (x[i] - mean_x) * (y[i] - mean_y);
    }

    return cov_sum / (double)(size - 1);
}

//Correlação de Pearson entre dois vetores

double cstats_correlation(const double *x, const double *y, size_t size) {
    if (!x || !y || size < 2) {
        return NAN;  // entrada inválida
    }

    // Calcular médias
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / (double)size;
    double mean_y = sum_y / (double)size;

    // Calcular covariância e variâncias
    double cov_sum = 0.0;
    double var_x_sum = 0.0;
    double var_y_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double dx = x[i] - mean_x;
        double dy = y[i] - mean_y;
        cov_sum += dx * dy;
        var_x_sum += dx * dx;
        var_y_sum += dy * dy;
    }

    double cov = cov_sum / (double)(size - 1);
    double stdev_x = sqrt(var_x_sum / (double)(size - 1));
    double stdev_y = sqrt(var_y_sum / (double)(size - 1));

    if (stdev_x == 0.0 || stdev_y == 0.0) {
        return NAN;  // evitar divisão por zero
    }

    return cov / (stdev_x * stdev_y);
}

/* =========================================================================
 * REGRESSÃO LINEAR
 * ========================================================================= */

//regressao linear simpples

 int cstats_linear_regression(const double *x, const double *y, size_t size,
                             double *slope, double *intercept) {
    if (!x || !y || !slope || !intercept || size < 2) {
        return -1;  // entrada inválida
    }

    // Calcular médias de x e y
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / (double)size;
    double mean_y = sum_y / (double)size;

    // Calcular soma dos produtos para o numerador e denominador
    double num = 0.0; // Σ (xi - mean_x)*(yi - mean_y)
    double den = 0.0; // Σ (xi - mean_x)^2
    for (size_t i = 0; i < size; i++) {
        double dx = x[i] - mean_x;
        num += dx * (y[i] - mean_y);
        den += dx * dx;
    }

    if (den == 0.0) {
        return -1;  // todos os x iguais → regressão indefinida
    }

    *slope = num / den;
    *intercept = mean_y - (*slope) * mean_x;

    return 0; // sucesso
}

//Regressão linear com métricas extras

int cstats_linear_regression_full(const double *x, const double *y, size_t size,
                                  double *slope,
                                  double *intercept,
                                  double *r_squared,
                                  double *residual_variance) {
    if (!x || !y || size < 2) {
        return -1;  // entrada inválida
    }

    // Calcular médias
    double sum_x = 0.0, sum_y = 0.0;
    for (size_t i = 0; i < size; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    double mean_x = sum_x / (double)size;
    double mean_y = sum_y / (double)size;

    // Calcular somas para slope
    double num = 0.0; // Σ (xi - mean_x)*(yi - mean_y)
    double den = 0.0; // Σ (xi - mean_x)^2
    for (size_t i = 0; i < size; i++) {
        double dx = x[i] - mean_x;
        num += dx * (y[i] - mean_y);
        den += dx * dx;
    }

    if (den == 0.0) {
        return -1;  // todos os x iguais → regressão indefinida
    }

    double a = num / den;              // slope
    double b = mean_y - a * mean_x;    // intercept

    if (slope) *slope = a;
    if (intercept) *intercept = b;

    // Calcular métricas extras se necessário
    double ss_total = 0.0;
    double ss_residual = 0.0;
    for (size_t i = 0; i < size; i++) {
        double y_pred = a * x[i] + b;
        double dy = y[i] - mean_y;
        double res = y[i] - y_pred;
        ss_total += dy * dy;
        ss_residual += res * res;
    }

    if (r_squared) *r_squared = (ss_total == 0.0) ? NAN : 1.0 - ss_residual / ss_total;
    if (residual_variance) *residual_variance = ss_residual / (double)(size - 2);

    return 0; // sucesso
}