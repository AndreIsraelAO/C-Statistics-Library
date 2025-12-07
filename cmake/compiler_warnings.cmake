set(COMMON_WARNINGS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wdouble-promotion
    -Wstrict-prototypes
)

add_compile_options(${COMMON_WARNINGS})
