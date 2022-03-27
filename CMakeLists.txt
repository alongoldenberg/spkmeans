cmake_minimum_required(VERSION 3.21)
project(spkmeans C)

set(CMAKE_C_STANDARD 11)


add_executable(spkmeans
        spkmeans.c)
target_link_libraries(spkmeans m)
target_link_libraries(spkmeans py)
target_link_libraries(spkmeans python3.9)
