cmake_minimum_required(VERSION 3.16)
project(Grid CXX)

set(CMAKE_CXX_STANDARD 17)

add_executable(simd_test simd_test.cpp
        SimdFieldTypes.h
        CMatOverloads.h
        LinkOverloads.h)

find_package(sleef CONFIG REQUIRED)
target_link_libraries(simd_test PRIVATE sleef::sleef)

target_compile_features(simd_test PRIVATE cxx_std_20)
target_compile_definitions(simd_test PRIVATE USE_SIMD_ACROSS_SITES)
