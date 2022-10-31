cmake_minimum_required(VERSION 3.16)
project(regulation_evolution)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -Wall)
set(CMAKE_CXX_FLAGS_DEBUG -g)
set(CMAKE_CXX_FLAGS_RELEASE -O2)

add_executable(cell_evolution
        src/cell_evolution.cpp
        src/ca.cpp
        src/hull.cpp
        src/gene.cpp
        src/genome.cpp
        src/cell.cpp
        src/conrec.cpp
        src/dish.cpp
        src/info.cpp
        src/misc.cpp
        src/output.cpp
        src/parameter.cpp
        src/parse.cpp
        src/pde.cpp
        src/intplane.cpp
        src/random.cpp
        src/crash.cpp
        src/warning.cpp
        src/x11graph.cpp
        src/ca.h
        src/hull.h
        src/gene.h
        src/genome.h
        src/cell.h
        src/conrec.h
        src/dish.h
        src/graph.h
        src/info.h
        src/misc.h
        src/output.h
        src/parameter.h
        src/parse.h
        src/pde.h
        src/intplane.h
        src/random.h
        src/sqr.h
        src/sticky.h
        src/crash.h
        src/warning.h
        src/x11graph.h src/boundingbox.cpp src/boundingbox.h)

find_package(Qt4 REQUIRED QtCore)
find_package(X11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

target_link_libraries(cell_evolution ${QT_LIBRARIES} X11 png)
target_compile_options(cell_evolution PRIVATE
        -Wno-unused-parameter
        -Wno-unused-variable
        -Wno-ignored-qualifiers
        -Wno-unused-result
        -Wno-write-strings
        -Wno-misleading-indentation)
