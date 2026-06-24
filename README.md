# C Statistics Library

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C Standard](https://img.shields.io/badge/C-C11-blue.svg)](https://en.wikipedia.org/wiki/C11_(C_standard_revision))
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

A robust, memory-safe, and numerically stable statistics library for C11. Provides opaque data structures for dynamic series management and a comprehensive suite of statistical functions with single-pass algorithms.

## Features

- **Memory-Safe Data Structures**: Opaque `StatSeries` type with automatic memory management and bounds checking
- **Numerical Stability**: Single-pass Welford's algorithm for variance and standard deviation
- **Comprehensive Statistics**: Mean, median, mode, quantiles, variance, correlation, and linear regression
- **Clean API**: Well-documented Doxygen-style headers for easy integration
- **Zero Dependencies**: Pure C11 with only standard library requirements (`math.h`, `stdlib.h`)

## Quick Start

Here's a minimal example demonstrating the core functionality:

```c
#include <stdio.h>
#include "stats.h"

int main(void) {
    // Create a new series with initial capacity of 10
    StatSeries* series = series_create(10);
    if (!series) {
        fprintf(stderr, "Failed to create series\n");
        return 1;
    }

    // Add values to the series
    series_add_value(series, 1.0f);
    series_add_value(series, 2.0f);
    series_add_value(series, 3.0f);
    series_add_value(series, 4.0f);
    series_add_value(series, 5.0f);

    // Access the data
    const float* data = series_data(series);
    size_t count = series_size(series);

    printf("Series contains %zu elements\n", count);
    printf("First element: %f\n", data[0]);

    // Clean up (always call destroy when done)
    series_destroy(series);

    return 0;
}
```

## Build & Installation

### Prerequisites

- GCC or Clang with C11 support
- GNU Make
- Standard C library with `math.h`

### Building from Source

```bash
# Clone the repository
git clone https://github.com/AndreIsraelAO/C-Statistics-Library.git
cd C-Statistics-Library

# Build the library and example
make clean && make all

# Run the example program
make example
```

### Installation

Install to system directories (requires sudo):

```bash
sudo make install
```

Custom installation prefix:

```bash
sudo make install PREFIX=/opt/custom
```

### Makefile Targets

| Target | Description |
|--------|-------------|
| `all` | Build static library (`libstats.a`) and example (default) |
| `example` | Build and run the example program |
| `install` | Install headers and library to PREFIX/include and PREFIX/lib |
| `clean` | Remove all build artifacts |
| `help` | Display available targets and variables |

## API Overview

### Dynamic Series Management

The `StatSeries` opaque type provides safe, dynamic array management:

```c
// Lifecycle
StatSeries* series_create(size_t initial_capacity);
void series_destroy(StatSeries* s);

// Data manipulation
bool series_add_value(StatSeries* s, float value);
size_t series_size(const StatSeries* s);
const float* series_data(const StatSeries* s);  // Read-only access
```

### Statistical Functions (cstats.h)

**Basic Metrics:**
- `cstats_mean()` - Arithmetic mean
- `cstats_median()` - Median value
- `cstats_mode()` / `cstats_multimode()` - Mode(s)

**Specialized Means:**
- `cstats_harmonic_mean()` - For rates and ratios
- `cstats_geometric_mean()` - For growth rates

**Dispersion:**
- `cstats_variance()` / `cstats_pvariance()` - Sample/population variance
- `cstats_stdev()` / `cstats_pstdev()` - Standard deviation

**Quantiles:**
- `cstats_quantile()` - Arbitrary quantile
- `cstats_quantiles()` - Multiple quantiles at once

**Relationships:**
- `cstats_covariance()` - Covariance between two series
- `cstats_correlation()` - Pearson correlation coefficient
- `cstats_linear_regression()` - Simple linear regression

## Project Structure

```
C-Statistics-Library/
├── include/           # Public headers
│   ├── stats.h        # StatSeries API
│   └── cstats.h       # Statistical functions
├── src/               # Implementation
│   ├── stats.c        # Dynamic series implementation
│   └── cstats.c       # Statistical algorithms
├── examples/          # Usage examples
│   └── example_basic.c
├── tests/             # Test suite
├── Makefile           # Build system
└── README.md          # This file
```

## Error Handling

All functions follow consistent error handling patterns:

- **Pointer functions**: Return `NULL` on allocation failure
- **Boolean functions**: Return `false` on error
- **Numeric functions**: Return `NaN` for invalid input (empty arrays, domain errors)
- **Integer functions**: Return `-1` or error code on failure

Always check return values before using results:

```c
StatSeries* s = series_create(10);
if (!s) {
    // Handle allocation failure
}

if (!series_add_value(s, value)) {
    // Handle reallocation failure
}
```

## Numerical Stability

The library implements **Welford's online algorithm** for computing variance in a single pass. This approach:

1. Avoids catastrophic cancellation from subtracting large similar numbers
2. Maintains numerical precision for large datasets
3. Processes data incrementally without storing all values

Formula update step for each new value `x`:
```
delta = x - mean_old
mean_new = mean_old + delta / n
M2_new = M2_old + delta * (x - mean_new)
variance = M2 / (n - 1)
```

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for details.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Support

- **Bug Reports**: [GitHub Issues](https://github.com/AndreIsraelAO/C-Statistics-Library/issues)
- **Discussions**: GitHub Discussions tab

---

Built with care for the C community. 📊
