cmake_minimum_required(VERSION 3.14)
project(tp1 C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(tp1
        neumann.c
        neumann.h
        skyline.c
        skyline.h
        testskyline.c)
