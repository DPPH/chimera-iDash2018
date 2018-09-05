cmake_minimum_required(VERSION 2.8)
project(fhe)

set(CMAKE_CXX_FLAGS "--std=c++17 -g3 -O0 -Wall -Werror -fopenmp")

set(CMAKE_CXX_FLAGS_DEBUG "--std=c++17 -g3 -O0 -Wall -Werror -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "--std=c++17 -g3 -O9 -DNDEBUG -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "-fopenmp")

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
        #BigFFT128.cpp
        #BigFFT128.h
        BigComplex.cpp
        BigComplex.h
        TLwe.cpp
        TLwe.h
        BigTorusVector.cpp BigTorusVector.h
        BigTorusPolynomial.cpp BigTorusPolynomial.h
        TRLwe.cpp TRLwe.h
        TRGSW.cpp TRGSW.h
        mainalgo.cpp mainalgo.h
        section2_params.h)

set(MAIN_PROGS
        main.cpp
        cloud-program.cpp
        cloud-program-debug.cpp
        key-gen-program.cpp
        debug-decrypt-program.cpp
        decrypt-program.cpp
        encrypt-program.cpp
        plaintext.cpp
        )

foreach (PROG_CPP ${MAIN_PROGS})
    get_filename_component(PROG ${PROG_CPP} NAME_WE)
    add_executable(${PROG} ${PROG_CPP})
    target_link_libraries(${PROG} chimera-dpph -lntl -lmpfr -lgmp -lntl -lm -pthread)
endforeach ()

add_executable(unittest
        test/arithmetic_test.cpp
        test/fft_test.cpp
        test/tlwe_test.cpp
        test/trlwe_test.cpp
        test/trgsw_test.cpp
        test/serialize_test.cpp
        test/mainalgo_test.cpp)
target_link_libraries(unittest gtest gtest_main chimera-dpph -lntl -lmpfr -lgmp -lm -pthread)

add_subdirectory(bench)

enable_testing()
add_test(unittest unittest)
