# Dont forget to set CMAKE_PREFIX_PATH to CONDA_PREFIX if dependencies live in conda environment
cmake_minimum_required(VERSION 3.16)
project(regulation_evolution)
include_directories(src)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -Wall)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(PNG REQUIRED)
find_package(X11 REQUIRED)

add_executable(cell_evolution
        src/cell_evolution.cpp
        src/ca.cpp src/ca.h
        src/hull.cpp src/hull.h
        src/gene.cpp src/gene.h
        src/genome.cpp src/genome.h
        src/cell.cpp src/cell.h
        src/conrec.cpp src/conrec.h
        src/dish.cpp src/dish.h
        src/misc.cpp src/misc.h
        src/output.cpp src/output.h
        src/parameter.cpp src/parameter.h
        src/parse.cpp src/parse.h
        src/pde.cpp src/pde.h
        src/intplane.cpp src/intplane.h
        src/random.cpp src/random.h
        src/crash.cpp
        src/warning.cpp src/warning.h
        src/x11graph.cpp src/x11graph.h
        src/graph.h
        src/sqr.h
        src/sticky.h
        src/boundingbox.cpp src/boundingbox.h
        src/foodpatch.cpp src/foodpatch.h)

target_link_libraries(cell_evolution PNG::PNG X11::X11)
target_compile_options(cell_evolution PRIVATE
        -Wno-unused-parameter
        -Wno-unused-variable
        -Wno-ignored-qualifiers
        -Wno-unused-result
        -Wno-write-strings
        -Wno-misleading-indentation)

add_executable(sweep_genomes
        scripts/genome/sweep_genomes.cpp
        src/genome.cpp src/genome.h
        src/gene.cpp src/gene.h
        src/misc.cpp src/misc.h
        src/random.cpp src/random.h)
