C-Statistics-Library

A C library for descriptive statistical calculations.
Provides simple, efficient, and portable functions for basic statistical metrics in C/C++ projects.

📌 Overview

C-Statistics-Library is a C11-compatible library that implements common descriptive statistics without external dependencies. It is designed to be easily integrated into:

Native C/C++ applications

Embedded systems

Mathematical backends exposed to other languages via FFI

The library focuses on simplicity, performance, and robustness, including validation for invalid inputs and NaN handling.

⚙️ Features

The library provides calculations for:

📊 Measures of Central Tendency

Arithmetic mean

Geometric mean

Harmonic mean

Median

Mode

📈 Measures of Dispersion

Variance

Standard deviation (population and sample)

Range

Quantiles

🔗 Correlation and Regression

Covariance

Pearson correlation coefficient

Linear regression (slope, intercept, and 
R2
R
2
)

Most functions run in O(n) time, except when sorting is required (e.g., median or mode).

📦 Installation
Requirements

C11-compatible compiler (GCC, Clang)

CMake 3.10+

Build instructions
git clone https://github.com/AndreIsraelAO/C-Statistics-Library.git
cd C-Statistics-Library

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make


Optional system-wide installation:

sudo make install

🚀 Quick Start
Basic example
#include <stdio.h>
#include <cstats.h>

int main() {
    double data[] = {10.5, 20.0, 45.2, 12.8, 10.5};
    size_t n = sizeof(data) / sizeof(data[0]);

    double mean = cstats_mean(data, n);
    double stdev = cstats_pstdev(data, n);

    printf("Size: %zu\n", n);
    printf("Mean: %.4f\n", mean);
    printf("Standard Deviation: %.4f\n", stdev);

    return 0;
}

Compile and run
gcc main.c -o stats_app -lcstats -lm
./stats_app

📘 Main API
Function	Description
cstats_mean	Arithmetic mean
cstats_median	Median
cstats_mode	Mode
cstats_variance	Variance
cstats_stdev	Standard deviation
cstats_covariance	Covariance
cstats_correlation	Pearson correlation
cstats_linear_regression	Linear regression

See the header file (include/cstats.h) for the complete API and function signatures.

🧪 Tests

Unit tests are available in the tests/ directory. To run them:

cd build
ctest

🤝 Contributing

Contributions are welcome. Suggested workflow:

Fork the repository

Create a feature branch:
git checkout -b feature/new-feature

Commit your changes

Open a Pull Request

Please follow the project’s coding style and ensure all tests pass.

📄 License

Distributed under the MIT License. See the LICENSE file for more information.

🧠 Recommendations for Improvement

To further professionalize the project:

Add usage examples or demos

Document all functions using Doxygen comments

Include CI badges (build status, test coverage)

Publish versioned releases with a clear changelog
