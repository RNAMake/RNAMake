cmake_minimum_required(VERSION 3.2)
project(RNAMake-Bindings)
set(CMAKE_CXX_STANDARD 17)
set(RNAMAKE /home/cjurich/projects/RNAMake/)

if(NOT DEFINED RNAMAKE)
    message(FATAL_ERROR "Variable RNAMAKE is not set. Please set it in CMakeLists.txt. Exiting...")
endif()

set(CMAKE_CXX_FLAGS -O3)

include_directories(src/)
include_directories(${RNAMAKE}/src/external/sqlite/)
add_compile_definitions(PYRNAMAKE_BINDINGS)

set(sql_src  ${RNAMAKE}/src/external/sqlite/sqlite3.c )

set_source_files_properties( ${RNAMAKE}/src/external/sqlite/sqlite3.h
		PROPERTIES
		LANGUAGE CXX )
set_source_files_properties( sql_src
		PROPERTIES
		LANGUAGE C
		LINKER CXX_STATIC_LIBRARY_LINKER )


add_library(sqlite3 ${RNAMAKE}/src/external/sqlite/sqlite3.c )
target_link_libraries(sqlite3 )


include_directories(${RNAMAKE}/src/)
include_directories(${RNAMAKE}/src/external/)

set(dirs
			base;
			data_structure;
			eternabot;
			math;
			motif;
			motif_data_structure;
			motif_search;
			motif_tools;
			resources;
			secondary_structure;
			sequence_optimization;
			structure;
			thermo_fluctuation;
			util;
			vienna
)
set(all_files)
foreach(dir ${dirs})
	file(GLOB_RECURSE temp ${RNAMAKE}/src/${dir}/*)
	list(APPEND all_files ${temp})
endforeach()

list(APPEND all_files ${sql_src})

# this is the python module
add_subdirectory(pybind11/)
# need to set the pybind11 directory
set(pybind11_DIR pybind11/)

file(GLOB_RECURSE module_includes src/*.hpp)

list(APPEND all_files ${module_includes})

add_library(all_lib ${all_files})
target_link_libraries( all_lib )

pybind11_add_module(RNAMake ${all_files} src/RNAMake.cpp)
target_link_libraries(RNAMake PUBLIC  all_lib sqlite3)
