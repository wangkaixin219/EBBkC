cmake_minimum_required(VERSION 3.6)
project(DESCol)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else ()
    # icpc
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
        message(icpc)
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        message(g++)
        set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    endif ()
endif ()


option(DUMP_COMPILATION_TO_FILE "dump compilation to files" ON)
if (NOT DUMP_COMPILATION_TO_FILE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
else ()
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)   # dump compilation commands to compile_commands.json
endif ()

add_compile_options(-O3 -g -march=native)
add_definitions(-DUSE_LOG -DTIME_RESULTS)


find_package(TBB REQUIRED)
include_directories(${TBB_INCLUDE_DIRS})
link_libraries(${TBB_LIBRARIES})

include_directories(${CMAKE_SOURCE_DIR}/truss/dependencies/sparsepp)
include_directories(${CMAKE_SOURCE_DIR}/truss/dependencies/libpopcnt)
include_directories(${CMAKE_SOURCE_DIR})


set(REORDER_UTILS truss/util/reordering/reorder_utils.cpp truss/util/reordering/pkc_reorder.cpp)
set(COMMON_UTILS truss/util/graph/graph.cpp truss/util/log/log.cpp truss/util/md5.cpp truss/util/search/search_util.cpp)
add_library(common-utils ${REORDER_UTILS} ${COMMON_UTILS} truss/decompose/extern_variables.cpp)
link_libraries(common-utils)
add_library(graph-pre-processing-legacy truss/util/graph/graph_preprocessing.cpp)
target_compile_definitions(graph-pre-processing-legacy PRIVATE NAVIE_SERIAL_CONSTRUCT=1)
add_library(graph-pre-processing truss/util/graph/graph_preprocessing.cpp)


#function(build_truss_algorithm name macros files)
#    add_executable(${name} ${files} truss/decompose/main.cpp truss/decompose/iter_helper.cpp truss/decompose/pkt_tri_cnt_intermediate.cpp)
#    target_link_libraries(${name} PRIVATE graph-pre-processing)
#
#    target_compile_definitions(${name} PRIVATE ${macros})
#    target_compile_definitions(${name} PRIVATE "BMP_PROCESSED=1")
#    target_compile_definitions(${name} PRIVATE "BMP_QUEUE=1")
#    #    target_compile_definitions(${name} PRIVATE "LOCAL_STAT=1") # switch off when do not collect local stat
#    #    target_compile_definitions(${name} PRIVATE "SUP_STAT=1") # switch off when do not collect local stat
#endfunction()

#build_truss_algorithm(truss "PKT_INTER=1" truss/decompose/pkt_merge.cpp)    ### Add BUCKET_BMP
#target_compile_definitions(truss PRIVATE "DISABLE_BUCKET_OPT=1;BUCKET_BMP=1")

add_executable(BBkC main.cpp def.cpp def.h edge_oriented.cpp edge_oriented.h set_operation.cpp set_operation.h
        truss/decompose/pkt_merge.cpp truss/decompose/iter_helper.cpp truss/decompose/pkt_tri_cnt_intermediate.cpp)
target_link_libraries(BBkC PRIVATE graph-pre-processing)
target_compile_definitions(BBkC PRIVATE "DISABLE_BUCKET_OPT=1;BUCKET_BMP=1;BMP_PROCESSED=1;BMP_QUEUE=1")

#target_compile_definitions(${name} PRIVATE ${macros})
#target_compile_definitions(${name} PRIVATE "BMP_PROCESSED=1")
#target_compile_definitions(${name} PRIVATE "BMP_QUEUE=1")







