cmake_minimum_required(VERSION 3.14)

project(
    cleanpt
    VERSION 0.1.0
    DESCRIPTION "A cleaner and modernized version of smallpt, a small path tracer"
    LANGUAGES CXX
)

# handle additional warnings
include(cmake/CompilerOptions.cmake)

# import sanitizers helpers
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)

add_executable(cleanpt main.cpp)
target_compile_features(cleanpt PRIVATE cxx_std_17)
add_warnings(cleanpt)
add_sanitizers(cleanpt)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(cleanpt PRIVATE OpenMP::OpenMP_CXX)
endif()

option(ENABLE_TESTING "Build tests" OFF)
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()