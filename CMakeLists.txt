cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(prog)

find_package(dlib REQUIRED)

add_executable(prog main.cpp)

target_link_libraries(prog PRIVATE dlib::dlib)
