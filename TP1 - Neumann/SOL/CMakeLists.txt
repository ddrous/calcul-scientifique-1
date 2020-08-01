cmake_minimum_required(VERSION 3.10)
project(project)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

file(GLOB project_SRC
    "*.h"
    "*.c"
)

add_executable(project ${project_SRC})


target_link_libraries(project m)