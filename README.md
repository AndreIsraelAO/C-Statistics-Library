# C-Statistics-Library

![Language](https://img.shields.io/badge/language-C11-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg) ![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey)

**CStatisticsLibrary** é uma biblioteca de estatística descritiva em C (C11), projetada para ser **simples, performática e segura**.

Focada em **determinismo e robustez**, ela trata `NaN`s e erros de entrada nativamente, sendo ideal para sistemas embarcados, kernels ou como backend matemático para linguagens de alto nível (Python, Rust, Go) via FFI.

---

## 🚀 Destaques

- **Zero Dependências:** Código C puro, sem bibliotecas externas.
- **Alta Performance:** Algoritmos otimizados com complexidade **O(n)** para a maioria das operações.
- **Safety First:** Tratamento robusto de erros (retorna `NAN` ou códigos de erro específicos).
- **FFI-Ready:** Estrutura de memória simples, facilitando *bindings* para outras linguagens.
- **Portável:** Compatível com C e C++, testado em GCC e Clang.

---

## 📦 Instalação

### Pré-requisitos
* Compilador C (GCC ou Clang)
* CMake 3.10+

### Compilando do código fonte (Universal)

```bash
git clone [https://github.com/AndreIsraelAO/C-Statistics-Library.git](https://github.com/AndreIsraelAO/C-Statistics-Library.git)
cd C-Statistics-Library
```

# Cria o diretório de build
```
mkdir build && cd build
```

# Configura e compila (Release para otimização máxima)
```
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

# Instalação (Opcional)
```
sudo make install
```

Para usuários Arch Linux

Se você (como eu) usa Arch, pode gerar e usar o pacote diretamente (adicione o PKGBUILD se tiver):
Bash

# Exemplo se você criar um PKGBUILD no futuro
```
makepkg -si
```

#⚡ Quick Start

Aqui está um exemplo completo de como usar a biblioteca no seu projeto.

```
#include <stdio.h>
#include <cstats.h> // Header principal

int main() {
    // Dataset de exemplo
    double dados[] = {10.5, 20.0, 45.2, 12.8, 10.5};
    size_t tamanho = sizeof(dados) / sizeof(dados[0]);

    // Cálculo de média
    double media = cstats_mean(dados, tamanho);
    
    // Cálculo de desvio padrão populacional
    double desvio = cstats_pstdev(dados, tamanho);

    printf("Dataset Size: %zu\n", tamanho);
    printf("Média: %.4f\n", media);
    printf("Desvio Padrão (Pop): %.4f\n", desvio);

    return 0;
}
```
Compilando o exemplo:
```

gcc main.c -o stats_app -lcstats -lm
./stats_app
```
#📚 Funcionalidades e API

A biblioteca cobre as principais métricas de estatística descritiva.
Categoria	Funções Disponíveis	Complexidade
Centralidade	mean, geometric_mean, harmonic_mean, median, mode	O(n) / O(n log n)*
Dispersão	variance, stdev, range, quantiles	O(n)
Correlação	covariance, correlation (Pearson)	O(n)
Regressão	linear_regression (Slope, Intercept, R²)	O(n)

*Mediana e modas podem exigir ordenação interna dependendo da implementação.

#🤝 Contribuindo

Contribuições são bem-vindas! Se você encontrou um bug ou quer melhorar a performance:

    Leia nosso Guia de Contribuição.

    Faça um fork do projeto.

    Crie sua Feature Branch (git checkout -b feature/NovaFeature).

    Commit suas mudanças (git commit -m 'Add some NovaFeature').

    Dê push para a branch (git push origin feature/NovaFeature).

    Abra um Pull Request.

#📄 Licença

Este projeto está licenciado sob a Licença MIT - veja o arquivo LICENSE para detalhes.

<p align="center"> Desenvolvido com 💙 e C por <a href="https://www.google.com/search?q=https://github.com/AndreIsraelAO">André Israel</a> </p>


