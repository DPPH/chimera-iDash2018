cmake_minimum_required(VERSION 2.8)
project(plaintext)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -g3 --std=c++11")

add_executable(plaintext main.cpp)
target_link_libraries(plaintext -lntl -lgmp -pthread)
