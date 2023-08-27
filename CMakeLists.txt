cmake_minimum_required(VERSION 3.0)

project(plain_elasticity)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable
add_executable(test main.cpp preprocess.cpp)

# Add include directories
target_include_directories(test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Add yaml-cpp library
find_library(YAML_CPP_LIB yaml-cpp PATHS /Users/chongran/lib/yaml-shared/lib REQUIRED)
target_include_directories(test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} /Users/chongran/lib/yaml-shared/include)
target_link_libraries(test PRIVATE ${YAML_CPP_LIB})