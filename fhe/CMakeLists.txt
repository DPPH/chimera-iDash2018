cmake_minimum_required(VERSION 3.10)
project(fhe)

set(CMAKE_CXX_STANDARD 17)

add_executable(fhe main.cpp BigFixP.cpp BigFixP.h)
target_link_libraries(fhe -lntl -lgmp -lm)
