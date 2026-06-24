#ifndef CSTATISTICS_STATISTICS_H
#define CSTATISTICS_STATISTICS_H

/**
 * @file cstats.h
 * @brief CStatisticsLibrary — Biblioteca de Estatística em C
 * 
 * Este arquivo define a API pública e estável da biblioteca.
 * As funções aqui declaradas são projetadas para serem:
 *   - Simples de usar
 *   - Confiáveis e determinísticas
 *   - Compatíveis com C e C++
 *   - Documentadas para facilitar manutenção futura
 */

#include <stddef.h>  /* para size_t */

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * MÉTRICAS BÁSICAS
 * ============================================================================ */

/**
 * @brief Calcula a média aritmética.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Média aritmética ou NaN se size == 0
 */
double cstats_mean(const double *data, size_t size);

/**
 * @brief Calcula a moda (valor mais frequente).
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Moda ou NaN se size == 0
 * 
 * @note Se houver múltiplas modas, retorna a primeira em ordem crescente.
 */
double cstats_mode(const double *data, size_t size);

/**
 * @brief Calcula todas as modas (multimode) de um array.
 * @param data       Ponteiro para os dados
 * @param size       Número de elementos
 * @param modes      Ponteiro para array que receberá as modas
 * @param max_modes  Tamanho máximo do array modes
 * @return Número de modas encontradas, ou -1 em caso de erro
 */
int cstats_multimode(const double *data, size_t size, 
                     double *modes, size_t max_modes);

/**
 * @brief Calcula a mediana.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Mediana ou NaN se size == 0
 */
double cstats_median(const double *data, size_t size);


/* ============================================================================
 * MÉDIAS VARIADAS
 * ============================================================================ */

/**
 * @brief Calcula a média harmônica.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Média harmônica ou NaN se algum valor for <= 0
 * 
 * Fórmula: n / (Σ 1/x_i)
 * @note Requer que todos os elementos sejam > 0.
 */
double cstats_harmonic_mean(const double *data, size_t size);

/**
 * @brief Calcula a média geométrica.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Média geométrica ou NaN se algum valor for <= 0
 * 
 * Fórmula: (Π x_i)^(1/n)
 * @note Requer que todos os elementos sejam > 0.
 */
double cstats_geometric_mean(const double *data, size_t size);


/* ============================================================================
 * MEDIANAS ESPECIALIZADAS
 * ============================================================================ */

/**
 * @brief Calcula a mediana baixa.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Mediana baixa ou NaN se size == 0
 * 
 * Em caso de número par de elementos, retorna o valor inferior do par central.
 * Ex.: [1, 2, 3, 4] → median_low = 2
 */
double cstats_median_low(const double *data, size_t size);

/**
 * @brief Calcula a mediana alta.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Mediana alta ou NaN se size == 0
 * 
 * Em caso de número par de elementos, retorna o valor superior do par central.
 * Ex.: [1, 2, 3, 4] → median_high = 3
 */
double cstats_median_high(const double *data, size_t size);

/**
 * @brief Calcula a mediana agrupada (grouped median).
 * @param data      Ponteiro para os dados
 * @param size      Número de elementos
 * @param interval  Tamanho do grupo (ex.: 1.0 para números inteiros)
 * @return Mediana agrupada ou NaN se size == 0
 * 
 * Fórmula: median = L + h * ((n/2 - cf) / f)
 * @note Segue o comportamento de statistics.median_grouped do Python.
 */
double cstats_median_grouped(const double *data, size_t size,
                             double interval);


/* ============================================================================
 * QUANTIS
 * ============================================================================ */

/**
 * @brief Calcula quantis arbitrários.
 * @param data      Ponteiro para os dados
 * @param size      Número de elementos
 * @param quantile  Valor entre 0.0 e 1.0
 * @return Quantil calculado ou NaN em caso de erro
 * 
 * Exemplos:
 *   - q = 0.5  → mediana
 *   - q = 0.25 → Q1 (primeiro quartil)
 *   - q = 0.75 → Q3 (terceiro quartil)
 */
double cstats_quantile(const double *data, size_t size, double quantile);

/**
 * @brief Calcula múltiplos quantis ao mesmo tempo.
 * @param data     Ponteiro para os dados
 * @param size     Número de elementos
 * @param q        Array de quantis desejados
 * @param q_count  Quantidade de quantis no array q
 * @param result   Array que receberá os resultados
 * @return 0 em sucesso, -1 em erro
 * 
 * Exemplo: q = {0.25, 0.5, 0.75}
 * @note result[] deve ter o mesmo tamanho de q_count.
 */
int cstats_quantiles(const double *data, size_t size,
                     const double *q, size_t q_count,
                     double *result);


/* ============================================================================
 * VARIÂNCIA E DESVIO PADRÃO
 * ============================================================================ */

/**
 * @brief Calcula a variância amostral.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Variância amostral ou NaN se size < 2
 * 
 * Fórmula: Σ (x - mean)² / (n - 1)
 */
double cstats_variance(const double *data, size_t size);

/**
 * @brief Calcula a variância populacional.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Variância populacional ou NaN se size < 1
 * 
 * Fórmula: Σ (x - mean)² / n
 */
double cstats_pvariance(const double *data, size_t size);

/**
 * @brief Calcula o desvio padrão amostral.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Desvio padrão amostral ou NaN se size < 2
 * 
 * Fórmula: sqrt(variance)
 */
double cstats_stdev(const double *data, size_t size);

/**
 * @brief Calcula o desvio padrão populacional.
 * @param data  Ponteiro para os dados
 * @param size  Número de elementos
 * @return Desvio padrão populacional ou NaN se size < 1
 * 
 * Fórmula: sqrt(pvariance)
 */
double cstats_pstdev(const double *data, size_t size);


/* ============================================================================
 * COVARIÂNCIA E CORRELAÇÃO
 * ============================================================================ */

/**
 * @brief Calcula a covariância amostral entre dois vetores.
 * @param x     Ponteiro para o primeiro vetor
 * @param y     Ponteiro para o segundo vetor
 * @param size  Número de elementos
 * @return Covariância amostral ou NaN se size < 2
 * 
 * Fórmula: Σ (xi - mean_x)(yi - mean_y) / (n - 1)
 */
double cstats_covariance(const double *x, const double *y, size_t size);

/**
 * @brief Calcula a correlação de Pearson entre dois vetores.
 * @param x     Ponteiro para o primeiro vetor
 * @param y     Ponteiro para o segundo vetor
 * @param size  Número de elementos
 * @return Coeficiente de correlação ou NaN se variâncias forem zero
 * 
 * Fórmula: cov(x, y) / (stdev(x) * stdev(y))
 */
double cstats_correlation(const double *x, const double *y, size_t size);


/* ============================================================================
 * REGRESSÃO LINEAR
 * ============================================================================ */

/**
 * @brief Realiza regressão linear simples (modelo: y = a*x + b).
 * @param x          Ponteiro para o vetor independente
 * @param y          Ponteiro para o vetor dependente
 * @param size       Número de elementos
 * @param slope      Ponteiro para receber o coeficiente angular (a)
 * @param intercept  Ponteiro para receber o intercepto (b)
 * @return 0 em sucesso, -1 em erro
 * 
 * @note Exige size >= 2.
 */
int cstats_linear_regression(const double *x, const double *y, size_t size,
                             double *slope, double *intercept);

/**
 * @brief Realiza regressão linear com métricas extras.
 * @param x                 Ponteiro para o vetor independente
 * @param y                 Ponteiro para o vetor dependente
 * @param size              Número de elementos
 * @param slope             Ponteiro para receber o coeficiente angular (opcional)
 * @param intercept         Ponteiro para receber o intercepto (opcional)
 * @param r_squared         Ponteiro para receber R² (opcional)
 * @param residual_variance Ponteiro para receber variância dos resíduos (opcional)
 * @return 0 em sucesso, -1 em erro
 * 
 * Calcula:
 *   - slope (coeficiente angular)
 *   - intercept (intercepto)
 *   - r² (coeficiente de determinação)
 *   - var_residuals (variância dos erros)
 * 
 * @note Qualquer ponteiro NULL é ignorado. Exige size >= 2.
 */
int cstats_linear_regression_full(const double *x, const double *y, size_t size,
                                  double *slope,
                                  double *intercept,
                                  double *r_squared,
                                  double *residual_variance);


#ifdef __cplusplus
}
#endif

#endif /* CSTATISTICS_STATISTICS_H */
