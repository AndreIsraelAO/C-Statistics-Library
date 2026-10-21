# C-Statistics-Library

[![CI](https://github.com/AndreIsraelAO/C-Statistics-Library/actions/workflows/ci.yml/badge.svg)](https://github.com/AndreIsraelAO/C-Statistics-Library/actions/workflows/ci.yml)
![License](https://img.shields.io/github/license/AndreIsraelAO/C-Statistics-Library)

**Languages:** English · [Português (Brasil)](README.pt-BR.md)

A small statistics library for C11 with no dependencies beyond the C standard library. The API is modeled on Python's `statistics` module, and the test suite checks results against reference values computed with it.

> **Status:** early-stage personal project (v0.1.0). It has unit tests and CI, but **no benchmarks** and no known production use. Performance has not been measured.

## Features

| Area | Functions |
| :--- | :--- |
| Central tendency | `mean`, `mode`, `multimode`, `median`, `median_low`, `median_high`, `median_grouped` |
| Other means | `harmonic_mean`, `geometric_mean`, `weighted_mean`, `trimmed_mean` |
| Quantiles | `quantile`, `quantiles` (linear interpolation, same as Python's `inclusive` method) |
| Dispersion | `variance`, `pvariance`, `stdev`, `pstdev`, `range`, `iqr`, `mad` |
| Shape | `skewness`, `kurtosis` (population moments; kurtosis is excess kurtosis) |
| Relationships | `covariance`, `correlation` (Pearson), `spearman` |
| Regression | `linear_regression`, `linear_regression_full` (slope, intercept, R², residual variance) |
| Distributions | Normal (`pdf`, `cdf`, `inv_cdf`), Binomial and Poisson (`pmf`, `cdf`) |
| Streaming | `cstats_accumulator`: Welford's algorithm for count, mean, variance, min and max, with `merge` to combine partial results; no heap allocation |

All functions are prefixed with `cstats_`. The full, documented API is in [`include/cstats.h`](include/cstats.h).

## Quick start

```c
#include <stdio.h>
#include <cstats.h>

int main(void) {
    double x[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[] = {2.1, 3.9, 6.2, 8.1, 10.5};
    double slope, intercept, r_squared, residual_var;

    printf("mean = %.2f, stdev = %.4f\n", cstats_mean(x, 5), cstats_stdev(x, 5));

    if (cstats_linear_regression_full(x, y, 5, &slope, &intercept,
                                      &r_squared, &residual_var) == 0) {
        printf("y = %.2fx + %.2f  (R² = %.4f)\n", slope, intercept, r_squared);
    }
    return 0;
}
```

Streaming statistics without storing the data:

```c
cstats_accumulator acc;
cstats_accumulator_init(&acc);
for (int i = 1; i <= 100; i++) cstats_accumulator_add(&acc, i);
printf("n=%zu mean=%.2f stdev=%.4f\n", cstats_accumulator_count(&acc),
       cstats_accumulator_mean(&acc), cstats_accumulator_stdev(&acc));
```

## Build and install

Requires a C11 compiler (GCC, Clang or MSVC) and CMake 3.16+.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build          # may need sudo
```

Options: `-DSTATS_BUILD_TESTS=ON|OFF` (default ON) and `-DSTATS_BUILD_SHARED=ON|OFF` (default OFF).

Using it from another CMake project after installing:

```cmake
find_package(cstats REQUIRED)
target_link_libraries(my_app PRIVATE cstats::cstats)
```

## Error handling

- Functions returning `double` return `NAN` on invalid input (`NULL` pointer, too few samples, parameter outside its domain) or when memory allocation fails.
- Functions returning `int` return `0` on success and `-1` on error.
- `cstats_multimode` returns the **total** number of modes (like `snprintf`); only the first `max_modes` are written.
- **NaN in the data:** functions that sort (median, mode, multimode, quantiles, grouped median, MAD, trimmed mean, Spearman) return `NAN`/`-1`, because ordering with NaN is undefined. Plain accumulations (mean, variance, covariance) propagate NaN per IEEE 754. `cstats_accumulator_add` rejects NaN without changing its state.
- `cstats_linear_regression_full` returns a `NAN` residual variance when `size == 2`.

## Design notes

- No global state: functions can be called from several threads on different data. This has not been verified with ThreadSanitizer.
- Functions that need ordering copy the input (`malloc`) and sort the copy with `qsort`, so the caller's array is never modified. They are O(n log n) and allocate. For streaming or allocation-free use, use the accumulator.
- Sums are plain `double` accumulation, without compensated (Kahan) summation.
- The public API is annotated with `CSTATS_API`, so the shared library exports symbols on Windows.

## Testing and CI

- `tests/test_main.c` is a dependency-free test runner (150 checks). Expected values come from Python's `statistics`, `math` and `NormalDist`.
- Edge cases covered: `NULL` and empty input, single samples, NaN in data, non-positive values in harmonic/geometric means, constant `x` in regression, `size == 2` residual variance, a too-small `multimode` buffer, out-of-domain distribution parameters, and agreement between the accumulator (including after `merge`) and the batch functions.
- GitHub Actions builds with warnings as errors (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wdouble-promotion`, or `/W4` on MSVC) on Linux (GCC, Clang), macOS (Clang) and Windows (MSVC), and runs one job under AddressSanitizer and UBSan.
- Locally verified with GCC (MinGW) and Clang on Windows; `find_package(cstats)` was checked with a separate consumer project.

## Known limitations and roadmap

- No hypothesis tests (t-test, chi-square, Kolmogorov–Smirnov) or confidence intervals yet; they need the incomplete gamma/beta functions and the Student's t distribution.
- Median and quantiles sort the whole array; a selection algorithm (quickselect) would make them O(n).
- Error causes are not distinguishable (`NAN` covers invalid input and allocation failure).
- Only `double` is supported.
- Planned: histograms and moving averages, `zscore`, compensated summation, benchmarks, Python bindings.

## Project layout

```text
include/cstats.h        public API and documentation
src/cstats.c            implementation
tests/test_main.c       test runner
cmake/                  warning flags helper
.github/workflows/      CI
```

## Contributing

Issues and pull requests are welcome. Please add a test with a reference value (ideally from Python's `statistics`) for any new function, and make sure the CI is green.

## License

MIT. See [LICENSE](LICENSE).

## Author

Andre Israel — [@AndreIsraelAO](https://github.com/AndreIsraelAO)
