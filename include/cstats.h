#ifndef CSTATS_H
#define CSTATS_H

/* Macro de exportação de símbolos (necessária para DLLs no Windows) */
#if defined(_WIN32) && !defined(CSTATS_STATIC)
#  ifdef CSTATS_BUILDING_DLL
#    define CSTATS_API __declspec(dllexport)
#  else
#    define CSTATS_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__)
#  define CSTATS_API __attribute__((visibility("default")))
#else
#  define CSTATS_API
#endif

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
CSTATS_API double cstats_mean(const double *data, size_t size);

/**
 * Calcula a moda (valor mais frequente).
 * Se houver múltiplas modas, retorna a primeira em ordem crescente.
 * Retorna NaN se size == 0.
 */
CSTATS_API double cstats_mode(const double *data, size_t size);

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
CSTATS_API int cstats_multimode(const double *data, size_t size, double *modes, size_t max_modes);


/**
 * Calcula a mediana.
 * Retorna NaN se size == 0.
 */
CSTATS_API double cstats_median(const double *data, size_t size);


/* =========================================================================
 * MÉDIAS VARIADAS
 * ========================================================================= */

/**
 * Média harmônica:
 *   n / (Σ 1/x_i)
 * Requer que todos os elementos sejam > 0.
 * Retorna NaN se algum valor for <= 0.
 */
CSTATS_API double cstats_harmonic_mean(const double *data, size_t size);

/**
 * Média geométrica:
 *   (Π x_i)^(1/n)
 * Requer que todos os elementos sejam > 0.
 * Retorna NaN se algum valor for <= 0.
 */
CSTATS_API double cstats_geometric_mean(const double *data, size_t size);


/* =========================================================================
 * MEDIANAS ESPECIALIZADAS
 * ========================================================================= */

/**
 * Mediana baixa:
 * Em caso de número par de elementos, retorna o valor inferior do par central.
 * Ex.: [1, 2, 3, 4] → median_low = 2
 */
CSTATS_API double cstats_median_low(const double *data, size_t size);

/**
 * Mediana alta:
 * Em caso de número par de elementos, retorna o valor superior do par central.
 * Ex.: [1, 2, 3, 4] → median_high = 3
 */
CSTATS_API double cstats_median_high(const double *data, size_t size);

/**
 * Mediana agrupada (grouped median).
 * Segue o comportamento de statistics.median_grouped do Python.
 *     median = L + h * ((n/2 - cf) / f)
 *
 * interval define o tamanho do grupo (ex.: 1.0 para números inteiros).
 * Retorna NaN em caso de size == 0.
 */
CSTATS_API double cstats_median_grouped(const double *data, size_t size,
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
CSTATS_API double cstats_quantile(const double *data, size_t size, double quantile);

/**
 * Atalho para múltiplos quantis ao mesmo tempo.
 * 
 * Exemplo: q = {0.25, 0.5, 0.75}
 * 
 * result[] deve ter o mesmo tamanho de q_count.
 * 
 * Retorna 0 em sucesso, -1 em erro.
 */
CSTATS_API int cstats_quantiles(const double *data, size_t size,
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
CSTATS_API double cstats_variance(const double *data, size_t size);

/**
 * Variância populacional.
 * Formula: Σ (x - mean)^2 / n
 * Requer n >= 1.
 */
CSTATS_API double cstats_pvariance(const double *data, size_t size);

/**
 * Desvio padrão amostral (sample standard deviation).
 * Igual a sqrt(variance).
 */
CSTATS_API double cstats_stdev(const double *data, size_t size);

/**
 * Desvio padrão populacional.
 * Igual a sqrt(pvariance).
 */
CSTATS_API double cstats_pstdev(const double *data, size_t size);


/* =========================================================================
 * COVARIÂNCIA E CORRELAÇÃO
 * ========================================================================= */

/**
 * Covariância amostral entre dois vetores.
 * Formula: Σ (xi - mean_x)*(yi - mean_y) / (n - 1)
 * Requer n >= 2.
 */
CSTATS_API double cstats_covariance(const double *x, const double *y, size_t size);

/**
 * Correlação de Pearson entre dois vetores.
 * Formula: cov(x, y) / (stdev(x) * stdev(y))
 * Retorna NaN se variâncias forem zero.
 */
CSTATS_API double cstats_correlation(const double *x, const double *y, size_t size);


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
CSTATS_API int cstats_linear_regression(const double *x, const double *y, size_t size,
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
CSTATS_API int cstats_linear_regression_full(const double *x, const double *y, size_t size,
                                  double *slope,
                                  double *intercept,
                                  double *r_squared,
                                  double *residual_variance);



/* =========================================================================
 * FORMA DA DISTRIBUIÇÃO
 * ========================================================================= */

/**
 * Assimetria (skewness) populacional de Fisher-Pearson:
 *   g1 = m3 / m2^(3/2), com m_k = Σ (x - mean)^k / n
 * Requer n >= 1. Retorna NaN se a variância for zero.
 */
CSTATS_API double cstats_skewness(const double *data, size_t size);

/**
 * Curtose em excesso (excess kurtosis) populacional:
 *   g2 = m4 / m2^2 - 3
 * Uma distribuição normal tem valor 0.
 * Requer n >= 1. Retorna NaN se a variância for zero.
 */
CSTATS_API double cstats_kurtosis(const double *data, size_t size);


/* =========================================================================
 * ACUMULADOR INCREMENTAL (STREAMING)
 * ========================================================================= */

/**
 * Acumulador online (algoritmo de Welford): calcula contagem, média,
 * variância, mínimo e máximo sem armazenar os dados e sem alocar memória.
 * Estável numericamente; dois acumuladores podem ser combinados com
 * cstats_accumulator_merge (útil para processamento paralelo).
 *
 * Os campos são privados: use apenas as funções abaixo.
 */
typedef struct cstats_accumulator {
    size_t count;
    double mean;
    double m2;   /* Σ (x - mean)^2 */
    double min;
    double max;
} cstats_accumulator;

/** Zera o acumulador. */
CSTATS_API void cstats_accumulator_init(cstats_accumulator *acc);

/** Adiciona uma observação. Retorna 0 em sucesso, -1 se acc for NULL ou x for NaN. */
CSTATS_API int cstats_accumulator_add(cstats_accumulator *acc, double x);

/** Combina src em dst (dst passa a representar a união das observações). */
CSTATS_API int cstats_accumulator_merge(cstats_accumulator *dst, const cstats_accumulator *src);

/** Número de observações (0 se acc for NULL). */
CSTATS_API size_t cstats_accumulator_count(const cstats_accumulator *acc);

/** Média. NaN se vazio. */
CSTATS_API double cstats_accumulator_mean(const cstats_accumulator *acc);

/** Variância amostral (n - 1). NaN se count < 2. */
CSTATS_API double cstats_accumulator_variance(const cstats_accumulator *acc);

/** Variância populacional (n). NaN se vazio. */
CSTATS_API double cstats_accumulator_pvariance(const cstats_accumulator *acc);

/** Desvio padrão amostral. NaN se count < 2. */
CSTATS_API double cstats_accumulator_stdev(const cstats_accumulator *acc);

/** Mínimo. NaN se vazio. */
CSTATS_API double cstats_accumulator_min(const cstats_accumulator *acc);

/** Máximo. NaN se vazio. */
CSTATS_API double cstats_accumulator_max(const cstats_accumulator *acc);


#ifdef __cplusplus
}
#endif

#endif /* CSTATS_H */
