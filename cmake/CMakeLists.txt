####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
project(RNAMake)
enable_testing()

include(/Users/jyesselm/projects/RNAMake/cmake/build/compiler.cmake)
set(RNAMAKE /Users/jyesselm/projects/RNAMake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories(/Users/jyesselm/projects/RNAMake/src/)
include_directories(/Users/jyesselm/projects/RNAMake/src/external/)
include_directories(/Users/jyesselm/projects/RNAMake/unittests/)
include_directories(/Users/jyesselm/projects/RNAMake/integration_tests)
include_directories(/Users/jyesselm/projects/RNAMake/apps/)
include(/Users/jyesselm/projects/RNAMake/cmake/build/sqlite.cmake)
 ####################################################################################################
####################################################################################################
# Library declarations 
####################################################################################################
###########################################################################
# base
###########################################################################
set(base_files
	/Users/jyesselm/projects/RNAMake/src/base/global_constants.cpp
	/Users/jyesselm/projects/RNAMake/src/base/log.cpp
	/Users/jyesselm/projects/RNAMake/src/base/exception.cpp
	/Users/jyesselm/projects/RNAMake/src/base/backtrace.cpp
	/Users/jyesselm/projects/RNAMake/src/base/types.hpp
	/Users/jyesselm/projects/RNAMake/src/base/paths.cpp
	/Users/jyesselm/projects/RNAMake/src/base/string.cpp 
)
add_library(base_lib  ${base_files})
target_link_libraries(base_lib  )
###########################################################################
# math
###########################################################################
set(math_files
	/Users/jyesselm/projects/RNAMake/src/math/rotation.cpp
	/Users/jyesselm/projects/RNAMake/src/math/vector_3.cpp
	/Users/jyesselm/projects/RNAMake/src/math/transform.cpp
	/Users/jyesselm/projects/RNAMake/src/math/stats.cpp
	/Users/jyesselm/projects/RNAMake/src/math/numerical.cpp
	/Users/jyesselm/projects/RNAMake/src/math/matrix_3x3.cpp 
)
add_library(math_lib  ${math_files})
target_link_libraries(math_lib base_lib  )
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	/Users/jyesselm/projects/RNAMake/src/data_structure/graph/graph.cc
	/Users/jyesselm/projects/RNAMake/src/data_structure/graph/graph_iter_list.cc
	/Users/jyesselm/projects/RNAMake/src/data_structure/graph/graph_base.cc
	/Users/jyesselm/projects/RNAMake/src/data_structure/graph/graph_adjacency_list.cc 
)
add_library(data_structure_lib  ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  )
###########################################################################
# util
###########################################################################
set(util_files
	/Users/jyesselm/projects/RNAMake/src/util/steric_lookup.cpp
	/Users/jyesselm/projects/RNAMake/src/util/x3dna/x3dna.cpp
	/Users/jyesselm/projects/RNAMake/src/util/find_pair.cpp
	/Users/jyesselm/projects/RNAMake/src/util/point_generator.cpp
)
add_library(util_lib  ${util_files})
target_link_libraries(util_lib math_lib sqlite3 )
add_library(all_lib  /Users/jyesselm/projects/RNAMake/src/main.cpp)
target_link_libraries(all_lib motif_tools_lib thermo_fluctuation_lib motif_search_lib sequence_optimization_lib )
####################################################################################################
####################################################################################################
# Unittest Declarations 
####################################################################################################
###########################################################################
# base Tests 
###########################################################################
	add_executable(string_unittest /Users/jyesselm/projects/RNAMake/unittests/base/string_unittest.cpp)
	target_link_libraries(string_unittest base_lib  )
	add_test(string_unittest string_unittest)
	add_executable(path_unittest /Users/jyesselm/projects/RNAMake/unittests/base/path_unittest.cpp)
	target_link_libraries(path_unittest base_lib  )
	add_test(path_unittest path_unittest)
	add_executable(log_unittest /Users/jyesselm/projects/RNAMake/unittests/base/log_unittest.cpp)
	target_link_libraries(log_unittest base_lib  )
	add_test(log_unittest log_unittest)
###########################################################################
# math Tests 
###########################################################################
	add_executable(vector_3_unittest /Users/jyesselm/projects/RNAMake/unittests/math/vector_3_unittest.cpp)
	target_link_libraries(vector_3_unittest math_lib  )
	add_test(vector_3_unittest vector_3_unittest)
	add_executable(transform_unittest /Users/jyesselm/projects/RNAMake/unittests/math/transform_unittest.cpp)
	target_link_libraries(transform_unittest math_lib  )
	add_test(transform_unittest transform_unittest)
	add_executable(rotation_unittest /Users/jyesselm/projects/RNAMake/unittests/math/rotation_unittest.cpp)
	target_link_libraries(rotation_unittest math_lib  )
	add_test(rotation_unittest rotation_unittest)
	add_executable(matrix_3x3_unittest /Users/jyesselm/projects/RNAMake/unittests/math/matrix_3x3_unittest.cpp)
	target_link_libraries(matrix_3x3_unittest math_lib  )
	add_test(matrix_3x3_unittest matrix_3x3_unittest)
###########################################################################
# data_structure Tests 
###########################################################################
	add_executable(graph_iter_list_unittest /Users/jyesselm/projects/RNAMake/unittests/data_structure/graph_iter_list_unittest.cpp)
	target_link_libraries(graph_iter_list_unittest data_structure_lib  )
	add_test(graph_iter_list_unittest graph_iter_list_unittest)
#	add_executable(graph_adjacency_list_unittest /Users/jyesselm/projects/RNAMake/unittests/data_structure/graph_adjacency_list_unittest.cpp)
#	target_link_libraries(graph_adjacency_list_unittest data_structure_lib  )
#	add_test(graph_adjacency_list_unittest graph_adjacency_list_unittest)
#	add_executable(graph_unittest /Users/jyesselm/projects/RNAMake/unittests/data_structure/graph_unittest.cpp)
#	target_link_libraries(graph_unittest data_structure_lib  )
#	add_test(graph_unittest graph_unittest)
#	add_executable(graph_base_unittest /Users/jyesselm/projects/RNAMake/unittests/data_structure/graph_base_unittest.cpp)
#	target_link_libraries(graph_base_unittest data_structure_lib  )
#	add_test(graph_base_unittest graph_base_unittest)
###########################################################################
# util Tests 
###########################################################################
#	add_executable(sqlite3_connection_unittest /Users/jyesselm/projects/RNAMake/unittests/util/sqlite3_connection_unittest.cpp)
#	target_link_libraries(sqlite3_connection_unittest util_lib  )
#	add_test(sqlite3_connection_unittest sqlite3_connection_unittest)
#	add_executable(motif_type_unittest /Users/jyesselm/projects/RNAMake/unittests/util/motif_type_unittest.cpp)
#	target_link_libraries(motif_type_unittest util_lib  )
#	add_test(motif_type_unittest motif_type_unittest)
	add_executable(steric_lookup_unittest /Users/jyesselm/projects/RNAMake/unittests/util/steric_lookup_unittest.cpp)
	target_link_libraries(steric_lookup_unittest util_lib  )
	add_test(steric_lookup_unittest steric_lookup_unittest)
	add_executable(find_pair_unittest /Users/jyesselm/projects/RNAMake/unittests/util/find_pair_unittest.cpp)
	target_link_libraries(find_pair_unittest util_lib  )
	add_test(find_pair_unittest find_pair_unittest)
#	add_executable(csv_unittest /Users/jyesselm/projects/RNAMake/unittests/util/csv_unittest.cpp)
#	target_link_libraries(csv_unittest util_lib  )
#	add_test(csv_unittest csv_unittest)
#	add_executable(uuid_unittest /Users/jyesselm/projects/RNAMake/unittests/util/uuid_unittest.cpp)
#	target_link_libraries(uuid_unittest util_lib  )
#	add_test(uuid_unittest uuid_unittest)
####################################################################################################
