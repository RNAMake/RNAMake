cmake_minimum_required(VERSION 3.2)
project(RNAMake-Bindings)
set(CMAKE_CXX_STANDARD 17)
set(RNAMAKE /home/cjurich/projects/RNAMake/)

if(NOT DEFINED RNAMAKE)
    message(FATAL_ERROR "Variable RNAMAKE is not set. Please set it in CMakeLists.txt. Exiting...")
endif()

include_directories(${RNAMAKE}/src/external/sqlite/)
add_compile_definitions(PYRNAMAKE_BINDINGS)

set(sql_src
		${RNAMAKE}/src/external/sqlite/sqlite3.c
		)

set_source_files_properties(
		${RNAMAKE}/src/external/sqlite/sqlite3.h
		PROPERTIES
		LANGUAGE CXX
)

set_source_files_properties(
		sql_src
		PROPERTIES
		LANGUAGE C
		LINKER CXX_STATIC_LIBRARY_LINKER
)


add_library(sqlite3 ${RNAMAKE}/src/external/sqlite/sqlite3.c )
target_link_libraries(sqlite3 )



include_directories(${RNAMAKE}/src/)
include_directories(${RNAMAKE}/src/external/)

set(dirs
base;
data_structure;
eternabot;
#external/sqlite;
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
	set(temp)
	file(GLOB_RECURSE temp ${RNAMAKE}/src/${dir}/*)
	foreach(file ${temp})
		list(APPEND all_files ${file})
	endforeach()

endforeach()

foreach(file ${sql_src})
	list(APPEND all_files ${file})
endforeach()

# this is the python module
add_subdirectory(pybind11/)
# need to set the pybind11 directory
set(pybind11_DIR pybind11/)

add_library(all_lib ${all_files})

include(bindings.cmake)
#target_link_libraries(all_lib)
#pybind11_add_module( bindings ${all_files} bindings.cpp )
#target_link_libraries( bindings PUBLIC all_lib )#sqlite3)
#pybind11_add_module( bindings2 ${all_files} bindings2.cpp )
#target_link_libraries( bindings2 PUBLIC all_lib )#sqlite3)
