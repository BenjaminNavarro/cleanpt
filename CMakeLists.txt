cmake_minimum_required(VERSION 3.14)

project(
    cleanpt
    VERSION 0.1.0
    DESCRIPTION "A cleaner and modernized version of smallpt, a small path tracer"
    LANGUAGES CXX
)

add_executable(cleanpt main.cpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cleanpt PRIVATE OpenMP::OpenMP_CXX)
endif()