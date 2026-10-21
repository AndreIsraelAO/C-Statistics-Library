#include "cstats.h"
#include <stddef.h>
#include <math.h>
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

// Cópia ordenada dos dados. Retorna NULL se houver NaN (a ordenação não
// seria bem definida) ou se a alocação falhar. O chamador deve usar free().
static double *sorted_copy(const double *data, size_t size) {
    double *copy = (double*)malloc(size * sizeof(double));
    if (!copy) return NULL;

    for (size_t i = 0; i < size; i++) {
        if (isnan(data[i])) {
            free(copy);
            return NULL;
        }
        copy[i] = data[i];
    }

    qsort(copy, size, sizeof(double), compare_doubles);
    return copy;
}

//mode
double cstats_mode(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;
    }

    double *sorted = sorted_copy(data, size);
    if (!sorted) {
        return NAN; // NaN nos dados ou falha de alocação
    }

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

    double *sorted = sorted_copy(data, size);
    if (!sorted) return -1; // NaN nos dados ou falha de alocação

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
                    modes[mode_index] = sorted[i - 1];
                }
                mode_index++;  // conta todas, mesmo as que não cabem em modes
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

    double *sorted = sorted_copy(data, size);
    if (!sorted) {
        return NAN;  // NaN nos dados ou falha de alocação
    }

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

    double *sorted = sorted_copy(data, size);
    if (!sorted) {
        return NAN; // NaN nos dados ou falha de alocação
    }

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

    double *sorted = sorted_copy(data, size);
    if (!sorted) {
        return NAN; // NaN nos dados ou falha de alocação
    }

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
    if (!data || size == 0 || !(interval > 0.0)) {
        return NAN;  // entrada inválida
    }

    double *sorted = sorted_copy(data, size);
    if (!sorted) return NAN;  // NaN nos dados ou falha de alocação

    // Mesmo algoritmo de statistics.median_grouped do Python
    double x = sorted[size / 2];
    size_t lo = 0;
    while (sorted[lo] < x) lo++;
    size_t hi = lo;
    while (hi < size && sorted[hi] == x) hi++;

    double cf = (double)lo;         // frequência acumulada antes de x
    double f = (double)(hi - lo);   // frequência de x
    double L = x - interval / 2.0;  // limite inferior do intervalo de x

    free(sorted);
    return L + interval * ((double)size / 2.0 - cf) / f;
}

/* =========================================================================
 * QUANTIS
 * ========================================================================= */

//quantil 

double cstats_quantile(const double *data, size_t size, double quantile) {
    if (!data || size == 0 || !(quantile >= 0.0 && quantile <= 1.0)) {
        return NAN;  // entrada inválida
    }

    double *sorted = sorted_copy(data, size);
    if (!sorted) return NAN;  // NaN nos dados ou falha de alocação

    double pos = quantile * (double)(size - 1);
    size_t index_lower = (size_t)floor(pos);
    size_t index_upper = (size_t)ceil(pos);
    double fraction = pos - (double)index_lower;

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

    double *sorted = sorted_copy(data, size);
    if (!sorted) return -1;  // NaN nos dados ou falha de alocação

    // Calcular cada quantil
    for (size_t i = 0; i < q_count; i++) {
        double qi = q[i];
        if (!(qi >= 0.0 && qi <= 1.0)) {
            free(sorted);
            return -1;  // quantil fora do intervalo
        }

        double pos = qi * (double)(size - 1);
        size_t index_lower = (size_t)floor(pos);
        size_t index_upper = (size_t)ceil(pos);
        double fraction = pos - (double)index_lower;

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

// Soma dos quadrados das diferenças em relação à média (size >= 1)
static double sum_sq_dev(const double *data, size_t size) {
    double mean = cstats_mean(data, size);

    double sq_diff_sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        sq_diff_sum += diff * diff;
    }
    return sq_diff_sum;
}

double cstats_variance(const double *data, size_t size) {
    if (!data || size < 2) {
        return NAN;  // entrada inválida
    }
    return sum_sq_dev(data, size) / (double)(size - 1);
}

//variancia populacional

double cstats_pvariance(const double *data, size_t size) {
    if (!data || size < 1) {
        return NAN;  // entrada inválida
    }
    return sum_sq_dev(data, size) / (double)size;
}

//desvio padrao amostral

double cstats_stdev(const double *data, size_t size) {
    return sqrt(cstats_variance(data, size));
}

//desvio padrao populacional

double cstats_pstdev(const double *data, size_t size) {
    return sqrt(cstats_pvariance(data, size));
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

    if (r_squared) *r_squared = (ss_total == 0.0) ? (double)NAN : 1.0 - ss_residual / ss_total;
    if (residual_variance) *residual_variance = (size > 2) ? ss_residual / (double)(size - 2) : (double)NAN;

    return 0; // sucesso
}

/* =========================================================================
 * FORMA DA DISTRIBUIÇÃO
 * ========================================================================= */

// Momento central de ordem k (k = 2, 3 ou 4), dividido por n
static double central_moment(const double *data, size_t size, double mean, int k) {
    double sum = 0.0;
    for (size_t i = 0; i < size; i++) {
        double diff = data[i] - mean;
        double term = diff * diff;
        if (k >= 3) term *= diff;
        if (k >= 4) term *= diff;
        sum += term;
    }
    return sum / (double)size;
}

double cstats_skewness(const double *data, size_t size) {
    if (!data || size < 1) {
        return NAN;  // entrada inválida
    }

    double mean = cstats_mean(data, size);
    double m2 = central_moment(data, size, mean, 2);
    if (m2 == 0.0) {
        return NAN;  // variância zero: assimetria indefinida
    }

    return central_moment(data, size, mean, 3) / pow(m2, 1.5);
}

double cstats_kurtosis(const double *data, size_t size) {
    if (!data || size < 1) {
        return NAN;  // entrada inválida
    }

    double mean = cstats_mean(data, size);
    double m2 = central_moment(data, size, mean, 2);
    if (m2 == 0.0) {
        return NAN;  // variância zero: curtose indefinida
    }

    return central_moment(data, size, mean, 4) / (m2 * m2) - 3.0;
}

/* =========================================================================
 * DISPERSÃO ROBUSTA E MÉDIAS PONDERADAS/APARADAS
 * ========================================================================= */

double cstats_range(const double *data, size_t size) {
    if (!data || size == 0) {
        return NAN;  // entrada inválida
    }

    double lo = data[0];
    double hi = data[0];
    for (size_t i = 0; i < size; i++) {
        if (isnan(data[i])) return NAN;
        if (data[i] < lo) lo = data[i];
        if (data[i] > hi) hi = data[i];
    }
    return hi - lo;
}

double cstats_iqr(const double *data, size_t size) {
    // NaN (entrada inválida ou NaN nos dados) se propaga pela subtração
    return cstats_quantile(data, size, 0.75) - cstats_quantile(data, size, 0.25);
}

double cstats_mad(const double *data, size_t size) {
    double med = cstats_median(data, size);
    if (isnan(med)) {
        return NAN;  // entrada inválida ou NaN nos dados
    }

    double *dev = (double*)malloc(size * sizeof(double));
    if (!dev) {
        return NAN;  // falha de alocação
    }
    for (size_t i = 0; i < size; i++) {
        dev[i] = fabs(data[i] - med);
    }

    double result = cstats_median(dev, size);
    free(dev);
    return result;
}

double cstats_weighted_mean(const double *data, const double *weights, size_t size) {
    if (!data || !weights || size == 0) {
        return NAN;  // entrada inválida
    }

    double weighted_sum = 0.0;
    double weight_total = 0.0;
    for (size_t i = 0; i < size; i++) {
        if (!(weights[i] >= 0.0)) {
            return NAN;  // pesos devem ser >= 0
        }
        weighted_sum += weights[i] * data[i];
        weight_total += weights[i];
    }

    if (!(weight_total > 0.0)) {
        return NAN;  // soma dos pesos deve ser positiva
    }
    return weighted_sum / weight_total;
}

double cstats_trimmed_mean(const double *data, size_t size, double proportion) {
    if (!data || size == 0 || !(proportion >= 0.0 && proportion < 0.5)) {
        return NAN;  // entrada inválida
    }

    double *sorted = sorted_copy(data, size);
    if (!sorted) {
        return NAN;  // NaN nos dados ou falha de alocação
    }

    size_t cut = (size_t)floor(proportion * (double)size);
    size_t kept = size - 2 * cut;  // sempre > 0 pois proportion < 0.5

    double sum = 0.0;
    for (size_t i = cut; i < cut + kept; i++) {
        sum += sorted[i];
    }

    free(sorted);
    return sum / (double)kept;
}

typedef struct {
    double value;
    size_t index;
} ranked_item;

static int compare_ranked(const void *a, const void *b) {
    double va = ((const ranked_item*)a)->value;
    double vb = ((const ranked_item*)b)->value;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

// Postos (1..n) com média dos postos em caso de empate. Retorna 0 ou -1.
static int average_ranks(const double *data, size_t size, double *ranks) {
    ranked_item *items = (ranked_item*)malloc(size * sizeof(ranked_item));
    if (!items) return -1;

    for (size_t i = 0; i < size; i++) {
        if (isnan(data[i])) {
            free(items);
            return -1;
        }
        items[i].value = data[i];
        items[i].index = i;
    }
    qsort(items, size, sizeof(ranked_item), compare_ranked);

    size_t i = 0;
    while (i < size) {
        size_t j = i;
        while (j + 1 < size && items[j + 1].value == items[i].value) j++;
        double rank = ((double)i + (double)j) / 2.0 + 1.0;
        for (size_t k = i; k <= j; k++) {
            ranks[items[k].index] = rank;
        }
        i = j + 1;
    }

    free(items);
    return 0;
}

double cstats_spearman(const double *x, const double *y, size_t size) {
    if (!x || !y || size < 2) {
        return NAN;  // entrada inválida
    }

    double *ranks = (double*)malloc(2 * size * sizeof(double));
    if (!ranks) {
        return NAN;  // falha de alocação
    }

    double result = NAN;
    if (average_ranks(x, size, ranks) == 0 &&
        average_ranks(y, size, ranks + size) == 0) {
        result = cstats_correlation(ranks, ranks + size, size);
    }

    free(ranks);
    return result;
}

/* =========================================================================
 * DISTRIBUIÇÕES DE PROBABILIDADE
 * ========================================================================= */

#define CSTATS_SQRT2 1.41421356237309504880
#define CSTATS_PI 3.14159265358979323846

double cstats_normal_pdf(double x, double mu, double sigma) {
    if (!(sigma > 0.0)) {
        return NAN;  // entrada inválida
    }

    double z = (x - mu) / sigma;
    return exp(-0.5 * z * z) / (sigma * sqrt(2.0 * CSTATS_PI));
}

double cstats_normal_cdf(double x, double mu, double sigma) {
    if (!(sigma > 0.0)) {
        return NAN;  // entrada inválida
    }

    double z = (x - mu) / sigma;
    return 0.5 * erfc(-z / CSTATS_SQRT2);
}

double cstats_normal_inv_cdf(double p, double mu, double sigma) {
    if (!(p > 0.0 && p < 1.0) || !(sigma > 0.0)) {
        return NAN;  // entrada inválida
    }

    // Bisseção sobre a CDF da normal padrão (monótona); converge até a
    // precisão de double em poucas iterações.
    double lo = -40.0;
    double hi = 40.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi) break;
        if (0.5 * erfc(-mid / CSTATS_SQRT2) < p) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return mu + sigma * (0.5 * (lo + hi));
}

double cstats_binomial_pmf(unsigned int k, unsigned int n, double p) {
    if (!(p >= 0.0 && p <= 1.0)) {
        return NAN;  // entrada inválida
    }
    if (k > n) return 0.0;
    if (p == 0.0) return k == 0 ? 1.0 : 0.0;
    if (p == 1.0) return k == n ? 1.0 : 0.0;

    double log_pmf = lgamma((double)n + 1.0) - lgamma((double)k + 1.0)
                   - lgamma((double)(n - k) + 1.0)
                   + (double)k * log(p) + (double)(n - k) * log1p(-p);
    return exp(log_pmf);
}

double cstats_binomial_cdf(unsigned int k, unsigned int n, double p) {
    if (!(p >= 0.0 && p <= 1.0)) {
        return NAN;  // entrada inválida
    }
    if (k >= n) return 1.0;

    double sum = 0.0;
    for (unsigned int i = 0; i <= k; i++) {
        sum += cstats_binomial_pmf(i, n, p);
    }
    return sum > 1.0 ? 1.0 : sum;
}

double cstats_poisson_pmf(unsigned int k, double lambda) {
    if (!(lambda >= 0.0) || isinf(lambda)) {
        return NAN;  // entrada inválida
    }
    if (lambda == 0.0) return k == 0 ? 1.0 : 0.0;

    return exp((double)k * log(lambda) - lambda - lgamma((double)k + 1.0));
}

double cstats_poisson_cdf(unsigned int k, double lambda) {
    if (!(lambda >= 0.0) || isinf(lambda)) {
        return NAN;  // entrada inválida
    }

    double sum = 0.0;
    for (unsigned int i = 0;; i++) {
        sum += cstats_poisson_pmf(i, lambda);
        if (i == k) break;
    }
    return sum > 1.0 ? 1.0 : sum;
}

/* =========================================================================
 * ACUMULADOR INCREMENTAL (STREAMING)
 * ========================================================================= */

void cstats_accumulator_init(cstats_accumulator *acc) {
    if (!acc) return;
    acc->count = 0;
    acc->mean = 0.0;
    acc->m2 = 0.0;
    acc->min = (double)NAN;
    acc->max = (double)NAN;
}

int cstats_accumulator_add(cstats_accumulator *acc, double x) {
    if (!acc || isnan(x)) {
        return -1;
    }

    acc->count++;
    double delta = x - acc->mean;
    acc->mean += delta / (double)acc->count;
    acc->m2 += delta * (x - acc->mean);

    if (acc->count == 1 || x < acc->min) acc->min = x;
    if (acc->count == 1 || x > acc->max) acc->max = x;
    return 0;
}

int cstats_accumulator_merge(cstats_accumulator *dst, const cstats_accumulator *src) {
    if (!dst || !src) {
        return -1;
    }
    if (src->count == 0) {
        return 0;
    }
    if (dst->count == 0) {
        *dst = *src;
        return 0;
    }

    // Combinação de Chan et al.
    double n_a = (double)dst->count;
    double n_b = (double)src->count;
    double n = n_a + n_b;
    double delta = src->mean - dst->mean;

    dst->m2 += src->m2 + delta * delta * n_a * n_b / n;
    dst->mean += delta * n_b / n;
    dst->count += src->count;
    if (src->min < dst->min) dst->min = src->min;
    if (src->max > dst->max) dst->max = src->max;
    return 0;
}

size_t cstats_accumulator_count(const cstats_accumulator *acc) {
    return acc ? acc->count : 0;
}

double cstats_accumulator_mean(const cstats_accumulator *acc) {
    if (!acc || acc->count == 0) return NAN;
    return acc->mean;
}

double cstats_accumulator_variance(const cstats_accumulator *acc) {
    if (!acc || acc->count < 2) return NAN;
    return acc->m2 / (double)(acc->count - 1);
}

double cstats_accumulator_pvariance(const cstats_accumulator *acc) {
    if (!acc || acc->count == 0) return NAN;
    return acc->m2 / (double)acc->count;
}

double cstats_accumulator_stdev(const cstats_accumulator *acc) {
    return sqrt(cstats_accumulator_variance(acc));
}

double cstats_accumulator_min(const cstats_accumulator *acc) {
    if (!acc || acc->count == 0) return NAN;
    return acc->min;
}

double cstats_accumulator_max(const cstats_accumulator *acc) {
    if (!acc || acc->count == 0) return NAN;
    return acc->max;
}