cmake_minimum_required(VERSION 3.16)

# Project metadata
project(
    CStatisticsLibrary
    VERSION 0.1.0
    LANGUAGES C
)

# C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Enable -fPIC (needed for shared libs, Python/Cython, etc.)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Options
option(STATS_BUILD_TESTS "Build unit tests" ON)
option(STATS_BUILD_SHARED "Build shared library" OFF)

# Compiler warnings (optional file)
include(cmake/compiler_warnings.cmake)

# Library type based on option
if(STATS_BUILD_SHARED)
    add_library(cstats SHARED src/cstats.c)
else()
    add_library(cstats STATIC src/cstats.c)
endif()

target_include_directories(cstats PUBLIC include)

########################
# Installation section #
########################

include(GNUInstallDirs)

install(
    TARGETS cstats
    EXPORT cstatsTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Export config for downstream use (C, Cython, CMake find_package)
install(
    EXPORT cstatsTargets
    FILE cstatsTargets.cmake
    NAMESPACE cstats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cstats
)

########################
# Testing section      #
########################

if(STATS_BUILD_TESTS)
    enable_testing()

    # Only the test source file(s) go into the executable — link with the cstats target
    add_executable(cstats_tests tests/test_main.c)
    target_link_libraries(cstats_tests PRIVATE cstats)

    add_test(NAME cstats_tests COMMAND cstats_tests)
endif()
