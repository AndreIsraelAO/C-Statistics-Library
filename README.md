# CStatisticsLibrary

**CStatisticsLibrary** é uma biblioteca de estatística em C, projetada para ser **simples, confiável e compatível com FFI** (Python, Rust, Go, etc).  

Ela fornece funções para **médias, medianas, moda, quantis, variância, desvio padrão, covariância, correlação e regressão linear**, com foco em **robustez, determinismo e performance**.

---

## Características principais

- API **pura em C**, fácil de compilar e usar.  
- Compatível com **C e C++**.  
- Projetada para **future bindings** (Python, Rust, Go).  
- Funções retornam `NAN` em entradas inválidas, garantindo **segurança e consistência**.  
- Implementação eficiente, com **complexidade O(n)** para a maioria das funções.  

---

## Funções Disponíveis
### Médias e Moda

double cstats_mean(const double *data, size_t size) – Média aritmética

double cstats_harmonic_mean(const double *data, size_t size) – Média harmônica

double cstats_geometric_mean(const double *data, size_t size) – Média geométrica

double cstats_mode(const double *data, size_t size) – Moda

int cstats_multimode(const double *data, size_t size, double *modes, size_t max_modes) – Todas as modas

### Medianas

double cstats_median(const double *data, size_t size) – Mediana

double cstats_median_low(const double *data, size_t size) – Mediana baixa

double cstats_median_high(const double *data, size_t size) – Mediana alta

double cstats_median_grouped(const double *data, size_t size, double interval) – Mediana agrupada

### Quantis

double cstats_quantile(const double *data, size_t size, double quantile) – Quantil único

int cstats_quantiles(const double *data, size_t size, const double *q, size_t q_count, double *result) – Quantis múltiplos

### Variância e Desvio Padrão

double cstats_variance(const double *data, size_t size) – Variância amostral

double cstats_pvariance(const double *data, size_t size) – Variância populacional

double cstats_stdev(const double *data, size_t size) – Desvio padrão amostral

double cstats_pstdev(const double *data, size_t size) – Desvio padrão populacional

### Covariância e Correlação

double cstats_covariance(const double *x, const double *y, size_t size) – Covariância amostral

double cstats_correlation(const double *x, const double *y, size_t size) – Correlação de Pearson

### Regressão Linear

int cstats_linear_regression(const double *x, const double *y, size_t size, double *slope, double *intercept) – Regressão linear simples

int cstats_linear_regression_full(const double *x, const double *y, size_t size, double *slope, double *intercept, double *r_squared, double *residual_variance) – Regressão linear com métricas completas

---

## Instalação

Clone o repositório e compile:

```bash
git clone https://github.com/seuusuario/C-Statistics-Library.git
cd C-Statistics-Library
mkdir build && cd build
cmake ..
make
sudo make install

