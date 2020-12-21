cmake_minimum_required(VERSION 3.2)
project(RNAMake-Bindings)
set(CMAKE_CXX_STANDARD 17)
set(RNAMAKE )

if(NOT DEFINED RNAMAKE)
    message(FATAL_ERROR "Variable RNAMAKE is not set. Please set it in CMakeLists.txt. Exiting...")
endif()

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

foreach(dir ${dirs})
	file(GLOB_RECURSE ${dir}_files ${RNAMAKE}/src/${dir}/*)
endforeach()

add_library(base_lib ${base_files})
target_link_libraries(base_lib)

add_library( math_lib ${math_files})
target_link_libraries(math_lib base_lib)
list(APPEND math_files ${base_files})

add_library(data_structure_lib ${data_structure_files})
target_link_libraries(data_structure_lib base_lib)
list(APPEND data_structure_files ${base_files})

add_library(util_lib ${util_files})
target_link_libraries(util_lib math_lib sqlite3)
list(APPEND util_files ${math_files})
list(APPEND util_files ${sql_src})

add_library(secondary_structure_lib ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib)
list(APPEND secondary_structure_files ${util_files})

add_library(eternabot_lib ${vienna_files} ${eternabot_files})
target_link_libraries(eternabot_lib secondary_structure_lib)
list(APPEND eternabot_files ${vienna_files})
list(APPEND eternabot_files ${secondary_structure_files})

add_library(structure_lib ${structure_files})
target_link_libraries(structure_lib util_lib)
list(APPEND structure_files ${util_files})

add_library(motif_lib ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib)
list(APPEND motif_files ${structure_files})
list(APPEND motif_files ${secondary_structure_files})

add_library(motif_tools_lib ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib)
list(APPEND motif_tools_files ${motif_files})

add_library(resources_lib ${resources_files})
target_link_libraries(resources_lib motif_lib)
list(APPEND resources_files ${motif_files})

add_library(motif_data_structure_lib ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib)
list(APPEND motif_data_structure_files ${resources_files})
list(APPEND motif_data_structure_files ${data_structure_files})

add_library(thermo_fluctuation_lib ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib)
list(APPEND thermo_fluctuation_files ${motif_data_structure_files})

add_library(motif_search_lib ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib)
list(APPEND motif_search_files ${motif_data_structure_files})

add_library(sequence_optimization_lib ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib eternabot_lib)
list(APPEND sequence_optimization_files ${motif_data_structure_files})
list(APPEND sequence_optimization_files ${eternabot_files})

# this is the python module
add_subdirectory(pybind11/)
# need to set the pybind11 directory
set(pybind11_DIR pybind11/)



pybind11_add_module( _base ${base_files} src/_base.cpp)
target_link_libraries( _base PUBLIC base_lib )

pybind11_add_module( _math ${math_files} src/_math.cpp)
target_link_libraries( _math PUBLIC math_lib )

pybind11_add_module( _data_structure ${data_structure_files} src/_data_structure.cpp)
target_link_libraries( _data_structure PUBLIC data_structure_lib )

pybind11_add_module( _util ${util_files} src/_util.cpp)
target_link_libraries( _util PUBLIC util_lib )

pybind11_add_module( _secondary_structure ${secondary_structure_files} src/_secondary_structure.cpp)
target_link_libraries( _secondary_structure PUBLIC secondary_structure_lib )

pybind11_add_module( _eternabot ${eternabot_files} src/_eternabot.cpp)
target_link_libraries( _eternabot PUBLIC eternabot_lib )

pybind11_add_module( _structure ${structure_files} src/_structure.cpp)
target_link_libraries( _structure PUBLIC structure_lib )

pybind11_add_module( _motif ${motif_files} src/_motif.cpp)
target_link_libraries( _motif PUBLIC motif_lib )

pybind11_add_module( _motif_tools ${motif_tools_files} src/_motif_tools.cpp)
target_link_libraries( _motif_tools PUBLIC motif_tools_lib )

pybind11_add_module( _resources ${resources_files} src/_resources.cpp)
target_link_libraries( _resources PUBLIC resources_lib )

pybind11_add_module( _motif_data_structure ${motif_data_structure_files} src/_motif_data_structure.cpp)
target_link_libraries( _motif_data_structure PUBLIC motif_data_structure_lib )


pybind11_add_module( _thermo_fluctuation ${thermo_fluctuation_files} src/_thermo_fluctuation.cpp)
target_link_libraries( _thermo_fluctuation PUBLIC thermo_fluctuation_lib )

pybind11_add_module( _motif_search ${motif_search_files} src/_motif_search.cpp)
target_link_libraries( _motif_search PUBLIC motif_search_lib )

pybind11_add_module( _sequence_optimization ${sequence_optimization_files} src/_sequence_optimization.cpp)
target_link_libraries( _sequence_optimization PUBLIC sequence_optimization_lib )
