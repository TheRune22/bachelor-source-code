cmake_minimum_required(VERSION 3.21)
project(bachelor-source-code)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)

add_compile_options(-Wall -Wextra -pedantic)

add_executable(main
        main.cpp algorithms.cpp algorithms.h graphs.cpp graphs.h utils.cpp utils.h experiments.cpp experiments.h)