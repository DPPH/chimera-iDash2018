cmake_minimum_required(VERSION 3.10)
project(plaintext)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Werror -g3")

set(CMAKE_CXX_STANDARD 17)

add_executable(plaintext main.cpp)
target_link_libraries(plaintext -lntl)
