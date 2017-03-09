cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
project (PhysModel)

# The version number.
set (PhysModel_VERSION_MAJOR 0)
set (PhysModel_VERSION_MINOR 0)

add_definitions(-std=c++11)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)
set (EXTRA_LIBS ${EXTRA_LIBS} src)

# add the executable
add_executable (PhysModel main.cpp)
target_link_libraries (PhysModel ${EXTRA_LIBS})
