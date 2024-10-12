cmake_minimum_required(VERSION 3.10)
project(ggq LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Boost debugging for more information
set(Boost_DEBUG 1)

# Find Boost with the required component
find_package(Boost 1.86 REQUIRED COMPONENTS unit_test_framework)

# Include directories for Boost
include_directories(${Boost_INCLUDE_DIRS})

# Add subdirectories (must be done after finding Boost)
add_subdirectory(src)
add_subdirectory(test)

# Add library
add_library(ggq_library 
    src/ggq.cpp 
    src/discretizer.cpp 
    src/utils/utils.cpp 
)

target_include_directories(ggq_library PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Link Boost::unit_test_framework with the library
target_link_libraries(ggq_library Boost::unit_test_framework)
