cmake_minimum_required(VERSION 2.8)
project(fhe)

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-g3 -O9 -DNDEBUG")

add_subdirectory(test/googletest)
include_directories(test/googletest/googletest)

add_library(chimera-dpph
        BigTorusFixP.cpp
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
        BigFFT128.cpp
        BigFFT128.h
        BigComplex.cpp
        BigComplex.h
        TLwe.cpp
        TLwe.h BigTorusVector.cpp BigTorusVector.h BigTorusPolynomial.cpp BigTorusPolynomial.h TRLwe.cpp TRLwe.h TRGSW.cpp TRGSW.h mainalgo.cpp mainalgo.h)

add_executable(fhe main.cpp )
target_link_libraries(fhe chimera-dpph -lntl -lmpfr -lgmp -lntl -lm -pthread)

add_executable(unittest test/arithmetic_test.cpp test/fft_test.cpp test/tlwe_test.cpp test/trlwe_test.cpp test/trgsw_test.cpp test/serialize_test.cpp test/mainalgo_test.cpp)
target_link_libraries(unittest gtest gtest_main chimera-dpph -lntl -lmpfr -lgmp -lm -pthread)

add_subdirectory(bench)

enable_testing()
add_test(unittest unittest)
