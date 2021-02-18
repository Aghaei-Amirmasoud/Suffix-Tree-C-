cmake_minimum_required(VERSION 3.16)
project(Suffixtree_0_2)

set(CMAKE_CXX_STANDARD 11)

add_executable(Suffixtree_0_2 main.cpp)

find_package(SFML 2.5.1 COMPONENTS system window graphics)

target_link_libraries(Suffixtree_0_2 sfml-system sfml-window sfml-graphics)


