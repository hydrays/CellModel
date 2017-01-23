cmake_minimum_required (VERSION 3.1.0 FATAL_ERROR)
project (PhysModel)

# The version number.
set (PhysModel_VERSION_MAJOR 0)
set (PhysModel_VERSION_MINOR 0)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories ("${PROJECT_BINARY_DIR}")

add_definitions(-std=c++11) 

include_directories ("${PROJECT_SOURCE_DIR}/src")
add_subdirectory (src)
set (EXTRA_LIBS ${EXTRA_LIBS} src)

#Check Boost
set(BOOST_MIN_VERSION "1.49.0")
set(BOOST_COMPONENTS date_time filesystem program_options system thread unit_test_framework random)
find_package(Boost ${BOOST_MIN_VERSION} COMPONENTS ${BOOST_COMPONENTS} REQUIRED)
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Fatal error: Boost (version >= 1.49.0) required.\n")
endif()
include_directories(${Boost_INCLUDE_DIR})

#find_path(GDAL "/usr/bin")

FIND_PACKAGE(GDAL REQUIRED)
include_directories(${GDAL_INCLUDE_DIR})
message(${GDAL_INCLUDE_DIR})
message(${GDAL_LIBRARIES})

#FIND_PACKAGE(Threads)
#FIND_PACKAGE(STXXL REQUIRED)

# add the executable
add_executable (PhysModel main.cpp)
target_link_libraries (PhysModel ${EXTRA_LIBS})
target_link_libraries (PhysModel ${Boost_LIBRARIES})
target_link_libraries (PhysModel ${GDAL_LIBRARIES})
target_link_libraries (PhysModel -lstxxl)
target_link_libraries (PhysModel -lpthread)
target_link_libraries (PhysModel -lcurl)
#target_link_libraries (PhysModel ${THREADS_LIBRARIES})

#target_compile_features(PhysModel PRIVATE cxx_range_for)
#set_property(TARGET PhysModel PROPERTY CXX_STANDARD 11)
#set_property(TARGET PhysModel PROPERTY CXX_STANDARD_REQUIRED ON)