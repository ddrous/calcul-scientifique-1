cmake_minimum_required(VERSION 3.14)
project(tp1_solution2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(tp1_solution2
        neumann.c
        neumann.h
        plotcom
        skyline.c
        skyline.h)

target_link_libraries(tp1_solution2 m)