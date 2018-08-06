cmake_minimum_required(VERSION 2.8)
project(fhe)

set(CMAKE_CXX_FLAGS "--std=c++17")

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
        BigComplex.h
        TLwe.cpp
        TLwe.h BigTorusVector.cpp BigTorusVector.h BigFixPVector.cpp BigFixPVector.h)

add_executable(fhe main.cpp )
target_link_libraries(fhe chimera-dpph -lntl -lmpfr -lgmp -lm -pthread)

add_executable(unittest test/arithmetic_test.cpp test/fft_test.cpp test/tlwe_test.cpp)
target_link_libraries(unittest gtest gtest_main chimera-dpph -lntl -lmpfr -lgmp -lm -pthread)

enable_testing()
add_test(unittest unittest)
