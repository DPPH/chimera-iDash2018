cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Werror")
set(CMAKE_CXX_FLAGS_RELEASE "-g3 -O9 -DNDEBUG")

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "No benchmark framework testing" FORCE)
add_subdirectory(benchmark EXCLUDE_FROM_ALL)

set(GOOGLEBENCH_SOURCES
        bench.cpp
        bench_fft.cpp
        bench_fft128.cpp
        )

#the benchmark applications are compiled with google benchmark framework
add_executable(bench ${GOOGLEBENCH_SOURCES})
target_link_libraries(bench benchmark chimera-dpph -lntl -lmpfr -lgmp -lm -pthread)
target_include_directories(bench
        PRIVATE benchmark/include
        PRIVATE ..
        )
