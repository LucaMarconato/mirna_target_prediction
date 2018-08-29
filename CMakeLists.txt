cmake_minimum_required (VERSION 3.6)
project (expression_profiles)
add_compile_options(-std=c++11 -Wall -g -O0)
# add_compile_options(-std=c++ -Ofast)
# add_definitions("-std=c++11 -Wall -Wextra -Wpedantic -g -O0")
# add_definitions("-std=c++11 -Ofast")

set (VERSION_MAJOR 0)
set (VERSION_MINOR 0)
set (VERSION_RELEASE 0)
set (VERSION_BUILD 0)

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY .)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE .)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG .)

include_directories(./dependencies /usr/local/Cellar/boost/1.67.0_1/include)
link_directories(/usr/local/Cellar/boost/1.67.0_1/lib)
file (GLOB SOURCES ./*.cpp)
add_executable(expression_profiles ${SOURCES})

add_custom_target(run COMMAND bash run.sh)
