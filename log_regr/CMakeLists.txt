cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-g3 -O9 -DNDEBUG")

set(TFHE_LIB_PATH external/tfhe_torus_generic/src)
set(CXXOPTS_LIB_PATH external/cxxopts)

add_subdirectory(${TFHE_LIB_PATH})
add_subdirectory(${CXXOPTS_LIB_PATH})

include_directories(
    ../common
    ${TFHE_LIB_PATH}/include
    ${CXXOPTS_LIB_PATH}/include
    )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(gen_keys gen_keys.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)

add_executable(encrypt_data encrypt_data.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)
target_link_libraries(encrypt_data -lntl -lgmp -pthread)

add_executable(decrypt_data decrypt_data.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)
target_link_libraries(decrypt_data -lntl -lgmp -pthread)

add_executable(log_regr log_regr.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)

add_executable(decrypt_beta decrypt_beta.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)
target_link_libraries(decrypt_beta -lntl -lgmp -pthread)

add_executable(decrypt_X_beta decrypt_X_beta.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)
target_link_libraries(decrypt_X_beta -lntl -lgmp -pthread)






add_subdirectory(external/googletest)
include_directories(external/googletest/googletest)

add_executable(unittest test/log_regr_fncs_test.cpp keyset.cpp
    $<TARGET_OBJECTS:tfhe-core>
    $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>
)
target_link_libraries(unittest gtest gtest_main -lntl -lgmp -pthread)

enable_testing()
add_test(unittest unittest)
