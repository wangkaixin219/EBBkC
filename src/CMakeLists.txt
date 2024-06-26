cmake_minimum_required(VERSION 3.16)
project(DESCol)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -mavx -fopenmp -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -mavx -fopenmp -O3")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
    message("ERROR: OpenMP could not be found.")
endif(OPENMP_FOUND)


add_executable(BBkC main.cpp
        def.cpp def.h
        edge_oriented.cpp edge_oriented.h
        set_operation.cpp set_operation.h)
