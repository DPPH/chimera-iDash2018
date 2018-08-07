cmake_minimum_required(VERSION 2.8)
project(plaintext)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -g3")

set(CMAKE_CXX_FLAGS "--std=c++17")


include_directories(../common)

add_executable(plaintext ${SRCS} main.cpp)
target_link_libraries(plaintext -lntl -lgmp -pthread)
