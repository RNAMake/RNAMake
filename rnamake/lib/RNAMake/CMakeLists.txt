cmake_minimum_required(VERSION 2.8.12)
project(rnamake_new)

set(CMAKE_BUILD_TYPE Release)
include(/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/cmake/build/compiler.cmake)

# Include path for Python header files
# Include paths for RNAMake src
include_directories(/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/)
include_directories(/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/plog/)
include_directories(/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/)

include_directories(/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/)

# sqlite libraries
find_library(SQLITE3_LIBRARY sqlite3)

###########################################################################
# base
###########################################################################

set(base_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/file_io.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/string.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/settings.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/application.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/cl_option.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/command_line_parser.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/global_constants.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/backtrace.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/log.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/base/option.cc
)

add_library(base_lib ${base_files})
target_link_libraries(base_lib)

add_executable(application_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/base/application_unittest.cpp)
target_link_libraries(application_unittest base_lib)

add_executable(option_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/base/option_unittest.cpp)
target_link_libraries(option_unittest base_lib)

add_executable(cl_option_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/base/cl_option_unittest.cpp)
target_link_libraries(cl_option_unittest base_lib)

add_executable(log_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/base/log_unittest.cc)
target_link_libraries(log_unittest base_lib)

###########################################################################
# math
###########################################################################

set(math_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/math/euler.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/math/hashing.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/math/stats.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/math/quaternion.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/math/numerical.cc
)

add_library(math_lib ${math_files})
target_link_libraries(math_lib base_lib )

add_executable(xyz_matrix_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/math/xyz_matrix_unittest.cpp)
target_link_libraries(xyz_matrix_unittest math_lib)

add_executable(quaternion_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/math/quaternion_unittest.cc)
target_link_libraries(quaternion_unittest math_lib)

add_executable(hashing_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/math/hashing_unittest.cc)
target_link_libraries(hashing_unittest math_lib)

###########################################################################
# data_structure
###########################################################################

set(data_structure_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph_adjacency_list.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph_base.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph_iter_list.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/tree/tree.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/tree/tree_node.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph/graph.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/data_structure/graph/graph_node.cc
)

add_library(data_structure_lib ${data_structure_files})
target_link_libraries(data_structure_lib base_lib )

add_executable(graph_adjacency_list_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/data_structure/graph_adjacency_list_unittest.cc)
target_link_libraries(graph_adjacency_list_unittest data_structure_lib)

add_executable(graph_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/data_structure/graph_unittest.cc)
target_link_libraries(graph_unittest data_structure_lib)

add_executable(tree_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/data_structure/tree_unittest.cpp)
target_link_libraries(tree_unittest data_structure_lib)

add_executable(graph_iter_list_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/data_structure/graph_iter_list_unittest.cc)
target_link_libraries(graph_iter_list_unittest data_structure_lib)

###########################################################################
# util
###########################################################################

set(util_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/cartesian_product.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/motif_type.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/uuid.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/sqlite3_connection.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/steric_lookup.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/monte_carlo.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/basic_io.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/util/x3dna.cc
)

add_library(util_lib ${util_files})
target_link_libraries(util_lib math_lib  ${SQLITE3_LIBRARY})

add_executable(uuid_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/util/uuid_unittest.cpp)
target_link_libraries(uuid_unittest util_lib)

add_executable(x3dna_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/util/x3dna_unittest.cpp)
target_link_libraries(x3dna_unittest util_lib)

add_executable(steric_lookup_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/util/steric_lookup_unittest.cpp)
target_link_libraries(steric_lookup_unittest util_lib)

add_executable(sqlite3_connection_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/util/sqlite3_connection_unittest.cpp)
target_link_libraries(sqlite3_connection_unittest util_lib)

###########################################################################
# vienna
###########################################################################

set(vienna_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/vienna/energy_par.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/vienna/vienna.cc
)

add_library(vienna_lib ${vienna_files})
target_link_libraries(vienna_lib base_lib )

###########################################################################
# secondary_structure
###########################################################################

set(secondary_structure_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/sequence_tools.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/util.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/secondary_structure_tree.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/structure.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/residue.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/motif.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/basepair.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/rna_structure.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/chain.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/pose.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/secondary_structure/secondary_structure_parser.cpp
)

add_library(secondary_structure_lib ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib )

add_executable(ss_pose_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/ss_pose_unittest.cpp)
target_link_libraries(ss_pose_unittest secondary_structure_lib)

add_executable(ss_structure_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/ss_structure_unittest.cpp)
target_link_libraries(ss_structure_unittest secondary_structure_lib)

add_executable(secondary_structure_parser_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/secondary_structure_parser_unittest.cpp)
target_link_libraries(secondary_structure_parser_unittest secondary_structure_lib)

add_executable(ss_residue_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/ss_residue_unittest.cpp)
target_link_libraries(ss_residue_unittest secondary_structure_lib)

add_executable(ss_motif_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/ss_motif_unittest.cpp)
target_link_libraries(ss_motif_unittest secondary_structure_lib)

add_executable(ss_chain_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/secondary_structure/ss_chain_unittest.cpp)
target_link_libraries(ss_chain_unittest secondary_structure_lib)

###########################################################################
# eternabot
###########################################################################

set(eternabot_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/sequence_designer.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/scorer.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/feature_generator.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy/num_of_yellow.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy/direction_of_gc.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy/berex_test.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy/clear_plot.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/eternabot/strategy/a_basic_test.cpp
)

add_library(eternabot_lib ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib  secondary_structure_lib )

###########################################################################
# structure
###########################################################################

set(structure_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/atom.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/residue_type_set.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/residue.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/basepair.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/structure.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/close_chain.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/pdb_parser.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/rna_structure.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/chain.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/is_equal.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/residue_type.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/structure/basepair_state.cc
)

add_library(structure_lib ${structure_files})
target_link_libraries(structure_lib util_lib )

add_executable(pdb_parser_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/pdb_parser_unittest.cpp)
target_link_libraries(pdb_parser_unittest structure_lib)

add_executable(structure_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/structure_unittest.cpp)
target_link_libraries(structure_unittest structure_lib)

add_executable(chain_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/chain_unittest.cpp)
target_link_libraries(chain_unittest structure_lib)

add_executable(residue_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/residue_unittest.cpp)
target_link_libraries(residue_unittest structure_lib)

add_executable(basepair_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/basepair_unittest.cpp)
target_link_libraries(basepair_unittest structure_lib)

add_executable(chain_closure_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/chain_closure_unittest.cc)
target_link_libraries(chain_closure_unittest structure_lib)

add_executable(atom_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/atom_unittest.cpp)
target_link_libraries(atom_unittest structure_lib)

add_executable(residue_type_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/residue_type_unittest.cpp)
target_link_libraries(residue_type_unittest structure_lib)

add_executable(basepair_state_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/structure/basepair_state_unittest.cpp)
target_link_libraries(basepair_state_unittest structure_lib)

###########################################################################
# motif
###########################################################################

set(motif_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_state_aligner.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_ensemble.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/pose.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_to_secondary_structure.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/pose_factory.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_scorer.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_state.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_state_ensemble.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif/motif_factory.cc
)

add_library(motif_lib ${motif_files})
target_link_libraries(motif_lib structure_lib  secondary_structure_lib )

add_executable(motif_ensemble_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_ensemble_unittest.cpp)
target_link_libraries(motif_ensemble_unittest motif_lib)

add_executable(motif_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_unittest.cpp)
target_link_libraries(motif_unittest motif_lib)

add_executable(motif_state_ensemble_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_state_ensemble_unittest.cpp)
target_link_libraries(motif_state_ensemble_unittest motif_lib)

add_executable(motif_factory_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_factory_unittest.cpp)
target_link_libraries(motif_factory_unittest motif_lib)

add_executable(motif_to_secondary_structure_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_to_secondary_structure_unittest.cpp)
target_link_libraries(motif_to_secondary_structure_unittest motif_lib)

add_executable(motif_state_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_state_unittest.cpp)
target_link_libraries(motif_state_unittest motif_lib)

add_executable(motif_scorer_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif/motif_scorer_unittest.cpp)
target_link_libraries(motif_scorer_unittest motif_lib)

###########################################################################
# motif_tools
###########################################################################

set(motif_tools_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_tools/segmenter.cpp
)

add_library(motif_tools_lib ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib )

add_executable(segmenter_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_tools/segmenter_unittest.cpp)
target_link_libraries(segmenter_unittest motif_tools_lib)

###########################################################################
# resources
###########################################################################

set(resources_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/motif_sqlite_library.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/motif_ensemble_sqlite_connection.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/sqlite_library.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/motif_state_ensemble_sqlite_library.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/added_motif_library.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/resource_manager.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/motif_sqlite_connection.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/resources/motif_state_sqlite_library.cpp
)

add_library(resources_lib ${resources_files})
target_link_libraries(resources_lib motif_lib )

add_executable(motif_sqlite_connection_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/motif_sqlite_connection_unittest.cpp)
target_link_libraries(motif_sqlite_connection_unittest resources_lib)

add_executable(resource_mananger_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/resource_mananger_unittest.cpp)
target_link_libraries(resource_mananger_unittest resources_lib)

add_executable(motif_sqlite_library_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/motif_sqlite_library_unittest.cpp)
target_link_libraries(motif_sqlite_library_unittest resources_lib)

add_executable(motif_ensemble_sqlite_connection_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/motif_ensemble_sqlite_connection_unittest.cpp)
target_link_libraries(motif_ensemble_sqlite_connection_unittest resources_lib)

add_executable(motif_state_sqlite_library_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/motif_state_sqlite_library_unittest.cpp)
target_link_libraries(motif_state_sqlite_library_unittest resources_lib)

add_executable(added_motif_library_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/added_motif_library_unittest.cpp)
target_link_libraries(added_motif_library_unittest resources_lib)

add_executable(motif_state_ensemble_sqlite_library_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/resources/motif_state_ensemble_sqlite_library_unittest.cpp)
target_link_libraries(motif_state_ensemble_sqlite_library_unittest resources_lib)

###########################################################################
# motif_data_structure
###########################################################################

set(motif_data_structure_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_state_tree.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_state_ensemble_tree.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_state_node.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_topology.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_merger.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_tree.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_state_graph.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_connection.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_graph.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_data_structure/motif_state_ensemble_graph.cc
)

add_library(motif_data_structure_lib ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib  data_structure_lib )

add_executable(motif_state_graph_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_state_graph_unittest.cpp)
target_link_libraries(motif_state_graph_unittest motif_data_structure_lib)

add_executable(motif_graph_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_graph_unittest.cpp)
target_link_libraries(motif_graph_unittest motif_data_structure_lib)

add_executable(motif_topology_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_topology_unittest.cpp)
target_link_libraries(motif_topology_unittest motif_data_structure_lib)

add_executable(motif_state_ensemble_tree_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_state_ensemble_tree_unittest.cpp)
target_link_libraries(motif_state_ensemble_tree_unittest motif_data_structure_lib)

add_executable(motif_merger_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_merger_unittest.cpp)
target_link_libraries(motif_merger_unittest motif_data_structure_lib)

add_executable(motif_state_tree_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_state_tree_unittest.cpp)
target_link_libraries(motif_state_tree_unittest motif_data_structure_lib)

add_executable(motif_tree_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_data_structure/motif_tree_unittest.cpp)
target_link_libraries(motif_tree_unittest motif_data_structure_lib)

###########################################################################
# thermo_fluctuation
###########################################################################

set(thermo_fluctuation_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/thermo_fluc_sampler.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/thermo_fluc_scorer.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/thermo_fluc_relax.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/graph/sampler.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/graph/scorer.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/thermo_fluctuation/graph/simulation.cc
)

add_library(thermo_fluctuation_lib ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib )

add_executable(thermo_fluc_sampler_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/thermo_fluctuation/thermo_fluc_sampler_unittest.cpp)
target_link_libraries(thermo_fluc_sampler_unittest thermo_fluctuation_lib)

add_executable(thermo_fluc_simulation_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/thermo_fluctuation/thermo_fluc_simulation_unittest.cpp)
target_link_libraries(thermo_fluc_simulation_unittest thermo_fluctuation_lib)

###########################################################################
# motif_search
###########################################################################

set(motif_search_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/motif_state_monte_carlo.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/problem.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/search.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/solution_topology.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/solution_filter.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/exhaustive/scorer.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/exhaustive/motif_state_enumerator.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/exhaustive/search.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/path_finding/node.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/path_finding/scorer.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/path_finding/search.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/path_finding/selector.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/monte_carlo/scorer.cc
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/motif_search/monte_carlo/search.cc
)

add_library(motif_search_lib ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib )

add_executable(monte_carlo_search_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_search/monte_carlo_search_unittest.cc)
target_link_libraries(monte_carlo_search_unittest motif_search_lib)

add_executable(exhaustive_search_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_search/exhaustive_search_unittest.cc)
target_link_libraries(exhaustive_search_unittest motif_search_lib)

add_executable(motif_state_search_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/motif_search/motif_state_search_unittest.cpp)
target_link_libraries(motif_state_search_unittest motif_search_lib)

###########################################################################
# sequence_optimization
###########################################################################

set(sequence_optimization_files
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/sequence_optimization/sequence_optimizer.cpp
	/Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/sequence_optimization/sequence_optimizer_3d.cpp
)

add_library(sequence_optimization_lib ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib  eternabot_lib )

add_executable(sequence_optimizer_unittest /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/unittests/sequence_optimization/sequence_optimizer_unittest.cpp)
target_link_libraries(sequence_optimizer_unittest sequence_optimization_lib)

###########################################################################
# all
###########################################################################

add_library(all_lib /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/src/main.cpp)
target_link_libraries(all_lib motif_tools_lib  thermo_fluctuation_lib  motif_search_lib  sequence_optimization_lib )

add_executable(sequence_optimizer_app /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/sequence_optimizer_app/sequence_optimizer_app.cpp )
target_link_libraries(sequence_optimizer_app all_lib)


add_executable(thermo_simulation /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/thermo_simulation/thermo_simulation.cpp )
target_link_libraries(thermo_simulation all_lib)


add_executable(sequence_optimization_benchmarks /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/sequence_optimization_benchmarks/sequence_optimization_benchmarks.cc )
target_link_libraries(sequence_optimization_benchmarks all_lib)


add_executable(design_rna_scaffold /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/design_rna_scaffold/design_rna_scaffold.cc )
target_link_libraries(design_rna_scaffold all_lib)


add_executable(hash_motifs /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/hash_motifs/hash_motifs.cc )
target_link_libraries(hash_motifs all_lib)


add_executable(get_best_solutions /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/get_best_solutions/get_best_solutions.cc )
target_link_libraries(get_best_solutions all_lib)


add_executable(eternabot /Users/jyesselm/projects/RNAMake/rnamake/lib/RNAMake/apps/eternabot/eternabot.cpp )
target_link_libraries(eternabot all_lib)


