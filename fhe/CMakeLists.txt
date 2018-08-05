cmake_minimum_required(VERSION 3.10)
project(fhe)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-g3 -O9 -DNDEBUG")

add_subdirectory(test/googletest)
include_directories(test/googletest/googletest)

add_library(chimera-dpph
        BigFixP.cpp
        BigFixP.h
        commons.cpp
        commons.h
        arithmetic.cpp
        arithmetic.h
        BigTorus.cpp
        BigTorus.h
        BigReal.cpp
        BigReal.h
        BigFFT.cpp
        BigFFT.h
        BigComplex.cpp
        BigComplex.h)

add_executable(fhe main.cpp )
target_link_libraries(fhe chimera-dpph -lntl -lmpfr -lgmp -lm)

add_executable(unittest test/arithmetic_test.cpp test/fft_test.cpp)
target_link_libraries(unittest gtest gtest_main chimera-dpph -lntl -lmpfr -lgmp -lm)

enable_testing()
add_test(unittest unittest)
