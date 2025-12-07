#ifndef CSTATISTICS_STATISTICS_H
#define CSTATISTICS_STATISTICS_H

/**
 * CStatisticsLibrary — Biblioteca de Estatística em C
 * ----------------------------------------------------
 * Este arquivo define a API pública e estável da biblioteca.
 * As funções aqui declaradas são projetadas para serem:
 *   - Simples de usar
 *   - Confiáveis e determinísticas
 *   - Compatíveis com C e C++
 *   - Documentadas para facilitar manutenção futura
 */

#include <stddef.h>  // para size_t

#ifdef __cplusplus
extern "C" {
#endif

/* =========================================================================
 * MÉTRICAS BÁSICAS
 * ========================================================================= */

/**
 * Calcula a média aritmética.
 * Retorna NaN se size == 0.
 */
double cstats_mean(const double *data, size_t size);

/**
 * Calcula a moda (valor mais frequente).
 * Se houver múltiplas modas, retorna a primeira em ordem crescente.
 * Retorna NaN se size == 0.
 */
double cstats_mode(const double *data, size_t size);

/**
 * Calcula todas as modas (multimode) de um array.
 *
 * param data      Ponteiro para os dados
 * param size      Número de elementos
 * param modes     Ponteiro para array que receberá as modas
 * param max_modes Tamanho máximo do array modes
 *
 * return Número de modas encontradas, ou -1 em caso de erro (size==0)
 */
int cstats_multimode(const double *data, size_t size, double *modes, size_t max_modes);


/**
 * Calcula a mediana.
 * Retorna NaN se size == 0.
 */
double cstats_median(const double *data, size_t size);


/* =========================================================================
 * MÉDIAS VARIADAS
 * ========================================================================= */

/**
 * Média harmônica:
 *   n / (Σ 1/x_i)
 * Requer que todos os elementos sejam > 0.
 * Retorna NaN se algum valor for <= 0.
 */
double cstats_harmonic_mean(const double *data, size_t size);

/**
 * Média geométrica:
 *   (Π x_i)^(1/n)
 * Requer que todos os elementos sejam > 0.
 * Retorna NaN se algum valor for <= 0.
 */
double cstats_geometric_mean(const double *data, size_t size);


/* =========================================================================
 * MEDIANAS ESPECIALIZADAS
 * ========================================================================= */

/**
 * Mediana baixa:
 * Em caso de número par de elementos, retorna o valor inferior do par central.
 * Ex.: [1, 2, 3, 4] → median_low = 2
 */
double cstats_median_low(const double *data, size_t size);

/**
 * Mediana alta:
 * Em caso de número par de elementos, retorna o valor superior do par central.
 * Ex.: [1, 2, 3, 4] → median_high = 3
 */
double cstats_median_high(const double *data, size_t size);

/**
 * Mediana agrupada (grouped median).
 * Segue o comportamento de statistics.median_grouped do Python.
 *     median = L + h * ((n/2 - cf) / f)
 *
 * interval define o tamanho do grupo (ex.: 1.0 para números inteiros).
 * Retorna NaN em caso de size == 0.
 */
double cstats_median_grouped(const double *data, size_t size,
                             double interval);

/* =========================================================================
 * QUANTIS
 * ========================================================================= */

/**
 * Calcula quantis arbitrários.
 * 
 * quantile ∈ [0.0, 1.0]
 * 
 * Exemplo:
 *   q=0.5 → mediana
 *   q=0.25 → Q1
 *   q=0.75 → Q3
 *
 * Retorna NaN se size == 0 ou se quantile estiver fora dos limites.
 */
double cstats_quantile(const double *data, size_t size, double quantile);

/**
 * Atalho para múltiplos quantis ao mesmo tempo.
 * 
 * Exemplo: q = {0.25, 0.5, 0.75}
 * 
 * result[] deve ter o mesmo tamanho de q_count.
 * 
 * Retorna 0 em sucesso, -1 em erro.
 */
int cstats_quantiles(const double *data, size_t size,
                     const double *q, size_t q_count,
                     double *result);


/* =========================================================================
 * VARIÂNCIA E DESVIO PADRÃO
 * ========================================================================= */

/**
 * Variância amostral (sample variance).
 * Formula: Σ (x - mean)^2 / (n - 1)
 * Requer n >= 2. Retorna NaN caso contrário.
 */
double cstats_variance(const double *data, size_t size);

/**
 * Variância populacional.
 * Formula: Σ (x - mean)^2 / n
 * Requer n >= 1.
 */
double cstats_pvariance(const double *data, size_t size);

/**
 * Desvio padrão amostral (sample standard deviation).
 * Igual a sqrt(variance).
 */
double cstats_stdev(const double *data, size_t size);

/**
 * Desvio padrão populacional.
 * Igual a sqrt(pvariance).
 */
double cstats_pstdev(const double *data, size_t size);


/* =========================================================================
 * COVARIÂNCIA E CORRELAÇÃO
 * ========================================================================= */

/**
 * Covariância amostral entre dois vetores.
 * Formula: Σ (xi - mean_x)*(yi - mean_y) / (n - 1)
 * Requer n >= 2.
 */
double cstats_covariance(const double *x, const double *y, size_t size);

/**
 * Correlação de Pearson entre dois vetores.
 * Formula: cov(x, y) / (stdev(x) * stdev(y))
 * Retorna NaN se variâncias forem zero.
 */
double cstats_correlation(const double *x, const double *y, size_t size);


/* =========================================================================
 * REGRESSÃO LINEAR
 * ========================================================================= */

/**
 * Regressão linear simples (modelo: y = a*x + b)
 *
 * Retorna:
 *   a → coeficiente angular (slope)
 *   b → intercepto (intercept)
 *
 * Exige size >= 2.
 * Retorna 0 em sucesso, -1 em erro.
 */
int cstats_linear_regression(const double *x, const double *y, size_t size,
                             double *slope, double *intercept);

/**
 * Regressão linear com métricas extras.
 *
 * Calcula:
 *   - slope
 *   - intercept
 *   - r²  (coeficiente de determinação)
 *   - var_residuals  (variância dos erros)
 *
 * Qualquer ponteiro NULL é ignorado.
 * size >= 2 é obrigatório.
 *
 * Retorna 0 em sucesso, -1 em erro.
 */
int cstats_linear_regression_full(const double *x, const double *y, size_t size,
                                  double *slope,
                                  double *intercept,
                                  double *r_squared,
                                  double *residual_variance);


#ifdef __cplusplus
}
#endif

#endif // CSTATISTICS_STATISTICS_H
