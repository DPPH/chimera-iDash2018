cmake_minimum_required(VERSION 3.12)
project(plaintext)

set(CMAKE_CXX_STANDARD 17)

add_executable(plaintext main.cpp)
target_link_libraries(plaintext -lntl)