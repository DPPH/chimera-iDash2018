cmake_minimum_required(VERSION 3.10)
project(fhe)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(test/googletest)
include_directories(test/googletest/googletest)

add_library(chimera-dpph
        BigFixP.cpp
        BigFixP.h
        commons.cpp
        commons.h arithmetic.cpp arithmetic.h BigTorus.cpp BigTorus.h)

add_executable(fhe main.cpp )
target_link_libraries(fhe chimera-dpph -lntl -lgmp -lm)

add_executable(unittest test/arithmetic_test.cpp)
target_link_libraries(unittest gtest gtest_main chimera-dpph -lntl -lgmp -lm)

enable_testing()
add_test(unittest unittest)
