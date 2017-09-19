cmake_minimum_required(VERSION 3.0)
project(PLAN)
set(PROJECT_NAME plan)

set(SOURCE_FILES  src/global.hpp src/global.cpp src/tree.hpp src/tree.cpp src/analyses.hpp src/analyses.cpp src/main.cpp)
add_executable(plan ${SOURCE_FILES})

cmake_host_system_information(RESULT myhostname QUERY HOSTNAME)
if (${myhostname} STREQUAL "Rixins-MacBook-Pro.local")
    # include MacPorts directories in my own laptop
    target_include_directories(${PROJECT_NAME} PUBLIC /opt/local/include /opt/local/include/openmpi-mp)
endif ()
############################################################
## Other users may add irregular header search paths here ##
############################################################

option(PARALLEL "Using parallel plan" OFF)
if (PARALLEL)
    message("============================================================")
    message("Generating Makefile for parallel program")
    message("============================================================")
    target_compile_options(${PROJECT_NAME} PUBLIC -std=c++11 -Wall -DMPI_ON -O3 -DNDEBUG -march=native)
    
    find_package(MPI REQUIRED)
    target_link_libraries(plan PUBLIC ${MPI_LIBRARIES})
else ()
    message("============================================================")
    message("Generating Makefile for serial program...")
    message("For parallel build, run cmake with -DPARALLEL=ON")
    message("============================================================")
    target_compile_options(${PROJECT_NAME} PUBLIC -std=c++11 -Wall -O3 -DNDEBUG -march=native)
endif ()

find_package(Boost REQUIRED)
target_link_libraries(plan PUBLIC ${Boost_LIBRARIES})

set(CMAKE_CXX_STANDARD 11)
