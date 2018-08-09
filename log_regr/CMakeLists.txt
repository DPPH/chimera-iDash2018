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

# set(LOG_REGR_SRCS
#         gen_keys.cpp
#         enc_log_regr.cpp
#     )

add_executable(gen_keys gen_keys.cpp $<TARGET_OBJECTS:tfhe-core> $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>)

add_executable(enc_log_regr enc_log_regr.cpp $<TARGET_OBJECTS:tfhe-core> $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>)
target_link_libraries(enc_log_regr -lntl -lgmp -pthread)
# target_include_directories(enc_log_regr
#         PRIVATE ${TFHE_LIB_PATH}/include
#         )


add_executable(log_regr log_regr.cpp $<TARGET_OBJECTS:tfhe-core> $<TARGET_OBJECTS:tfhe-fft-spqlios-avx>)
