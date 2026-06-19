cmake_minimum_required(VERSION 3.14)
project(PulseForge LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra)

file(GLOB_RECURSE BSP_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Core DSP / biosignal library.
add_library(bsp STATIC ${BSP_SOURCES})
target_include_directories(bsp PUBLIC ${CMAKE_SOURCE_DIR}/include)

# Command-line application.
add_executable(pulseforge app/main.cpp)
target_link_libraries(pulseforge PRIVATE bsp)

# Tests.
enable_testing()
add_executable(bsp_tests tests/test_main.cpp)
target_link_libraries(bsp_tests PRIVATE bsp)
add_test(NAME bsp_tests COMMAND bsp_tests)
