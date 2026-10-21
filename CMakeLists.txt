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

# Compiler warnings (per-target helper)
include(cmake/compiler_warnings.cmake)

# Library type based on option
if(STATS_BUILD_SHARED)
    add_library(cstats SHARED src/cstats.c)
    target_compile_definitions(cstats PRIVATE CSTATS_BUILDING_DLL)
else()
    add_library(cstats STATIC src/cstats.c)
    target_compile_definitions(cstats PUBLIC CSTATS_STATIC)
endif()

cstats_enable_warnings(cstats)

target_include_directories(cstats
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# libm is required for log/exp/sqrt/floor/ceil on UNIX-like systems (not on Windows)
if(UNIX)
    target_link_libraries(cstats PUBLIC m)
endif()

########################
# Installation section #
########################

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

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

# Export config for downstream use (find_package(cstats))
install(
    EXPORT cstatsTargets
    FILE cstatsTargets.cmake
    NAMESPACE cstats::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cstats
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cstatsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cstatsConfig.cmake
"include(\"\${CMAKE_CURRENT_LIST_DIR}/cstatsTargets.cmake\")\n")

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cstatsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cstatsConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cstats
)

########################
# Testing section      #
########################

if(STATS_BUILD_TESTS)
    enable_testing()

    add_executable(cstats_tests tests/test_main.c)
    target_link_libraries(cstats_tests PRIVATE cstats)
    cstats_enable_warnings(cstats_tests)

    add_test(NAME cstats_tests COMMAND cstats_tests)
endif()
