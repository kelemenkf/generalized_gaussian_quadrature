cmake_minimum_required(VERSION 3.10)
project(ggq LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(test)

add_library(ggq_library src/ggq.cpp src/discretizer.cpp src/utils/utils.cpp include/ggq.hpp include/discretizer.hpp include/utils.hpp)

