####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
project(RNAMake)

include(/Users/cjurich/projects/RNAMake/cmake/build/compiler.cmake)
set(RNAMAKE /Users/cjurich/projects/RNAMake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories(/Users/cjurich/projects/RNAMake/src/)
include_directories(/Users/cjurich/projects/RNAMake/src/external/)
include_directories(/Users/cjurich/projects/RNAMake/unittests/)
include_directories(/Users/cjurich/projects/RNAMake/apps/)
 ####################################################################################################
####################################################################################################
# Library declarations 
####################################################################################################
###########################################################################
# base
###########################################################################
set(base_files
	/Users/cjurich/projects/RNAMake/src/base/backtrace.cpp
	/Users/cjurich/projects/RNAMake/src/base/env_manager.cpp
	/Users/cjurich/projects/RNAMake/src/base/string.cc
	/Users/cjurich/projects/RNAMake/src/base/sys_interface.cpp
	/Users/cjurich/projects/RNAMake/src/base/vector_container.h
	/Users/cjurich/projects/RNAMake/src/base/command_line_parser.cpp
	/Users/cjurich/projects/RNAMake/src/base/application.cpp
	/Users/cjurich/projects/RNAMake/src/base/log.cc
	/Users/cjurich/projects/RNAMake/src/base/cl_option.cc
	/Users/cjurich/projects/RNAMake/src/base/file_io.cc
	/Users/cjurich/projects/RNAMake/src/base/assertions.h
	/Users/cjurich/projects/RNAMake/src/base/global_constants.cc
	/Users/cjurich/projects/RNAMake/src/base/util.hpp
	/Users/cjurich/projects/RNAMake/src/base/settings.cc
	/Users/cjurich/projects/RNAMake/src/base/types.h
	/Users/cjurich/projects/RNAMake/src/base/option.cc 
)
add_library(base_lib  ${base_files})
target_link_libraries(base_lib  )
###########################################################################
# math
###########################################################################
set(math_files
	/Users/cjurich/projects/RNAMake/src/math/euler.cc
	/Users/cjurich/projects/RNAMake/src/math/quaternion.cc
	/Users/cjurich/projects/RNAMake/src/math/xyz_vector.h
	/Users/cjurich/projects/RNAMake/src/math/xyz_matrix.h
	/Users/cjurich/projects/RNAMake/src/math/transform.h
	/Users/cjurich/projects/RNAMake/src/math/stats.cc
	/Users/cjurich/projects/RNAMake/src/math/hashing.cc
	/Users/cjurich/projects/RNAMake/src/math/numerical.cc 
)
add_library(math_lib  ${math_files})
target_link_libraries(math_lib base_lib  )
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	/Users/cjurich/projects/RNAMake/src/data_structure/graph_adjacency_list.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph/graph.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph_base.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph_iter_list.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph/graph_node.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/graph/graph_node.fwd.h
	/Users/cjurich/projects/RNAMake/src/data_structure/tree/tree_node.cc
	/Users/cjurich/projects/RNAMake/src/data_structure/tree/tree.cc 
)
add_library(data_structure_lib  ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  )
###########################################################################
# util
###########################################################################
set(util_files
	/Users/cjurich/projects/RNAMake/src/util/random_number_generator.h
	/Users/cjurich/projects/RNAMake/src/util/dssr.cpp
	/Users/cjurich/projects/RNAMake/src/util/steric_lookup.cpp
	/Users/cjurich/projects/RNAMake/src/util/basic_io.cpp
	/Users/cjurich/projects/RNAMake/src/util/cartesian_product.cc
	/Users/cjurich/projects/RNAMake/src/util/csv.cpp
	/Users/cjurich/projects/RNAMake/src/util/x3dna.cc
	/Users/cjurich/projects/RNAMake/src/util/monte_carlo.cpp
	/Users/cjurich/projects/RNAMake/src/util/motif_type.cc
	/Users/cjurich/projects/RNAMake/src/util/sqlite3_connection.cc
	/Users/cjurich/projects/RNAMake/src/util/uuid.cc 
)
add_library(util_lib  ${util_files})
target_link_libraries(util_lib math_lib sqlite3 )
###########################################################################
# vienna
###########################################################################
set(vienna_files
	/Users/cjurich/projects/RNAMake/src/vienna/pair_mat.h
	/Users/cjurich/projects/RNAMake/src/vienna/energy_par.cc
	/Users/cjurich/projects/RNAMake/src/vienna/vienna.cc 
)
add_library(vienna_lib  ${vienna_files})
target_link_libraries(vienna_lib base_lib  )
###########################################################################
# secondary_structure
###########################################################################
set(secondary_structure_files
	/Users/cjurich/projects/RNAMake/src/secondary_structure/util.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/sequence_constraint.cc
	/Users/cjurich/projects/RNAMake/src/secondary_structure/basepair.cc
	/Users/cjurich/projects/RNAMake/src/secondary_structure/rna_structure.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/secondary_structure_parser.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/chain.cc
	/Users/cjurich/projects/RNAMake/src/secondary_structure/sequence_tools.cc
	/Users/cjurich/projects/RNAMake/src/secondary_structure/structure.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/secondary_structure_tree.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/motif.cc
	/Users/cjurich/projects/RNAMake/src/secondary_structure/pose.cpp
	/Users/cjurich/projects/RNAMake/src/secondary_structure/residue.cc 
)
add_library(secondary_structure_lib  ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib  )
###########################################################################
# eternabot
###########################################################################
set(eternabot_files
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy/berex_test.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy/direction_of_gc.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/scorer.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy/clear_plot.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy/a_basic_test.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/feature_generator.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/sequence_designer.cpp
	/Users/cjurich/projects/RNAMake/src/eternabot/strategy/num_of_yellow.cpp 
)
add_library(eternabot_lib  ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib secondary_structure_lib  )
###########################################################################
# structure
###########################################################################
set(structure_files
	/Users/cjurich/projects/RNAMake/src/structure/is_equal.cpp
	/Users/cjurich/projects/RNAMake/src/structure/pdb_parser.cc
	/Users/cjurich/projects/RNAMake/src/structure/residue_type.cc
	/Users/cjurich/projects/RNAMake/src/structure/cif_parser.cc
	/Users/cjurich/projects/RNAMake/src/structure/structure.cc
	/Users/cjurich/projects/RNAMake/src/structure/residue.cc
	/Users/cjurich/projects/RNAMake/src/structure/chain.fwd.h
	/Users/cjurich/projects/RNAMake/src/structure/basepair.cc
	/Users/cjurich/projects/RNAMake/src/structure/exceptions.cc
	/Users/cjurich/projects/RNAMake/src/structure/close_chain.cc
	/Users/cjurich/projects/RNAMake/src/structure/atom.cc
	/Users/cjurich/projects/RNAMake/src/structure/chain.cc
	/Users/cjurich/projects/RNAMake/src/structure/basepair_state.cc
	/Users/cjurich/projects/RNAMake/src/structure/rna_structure.cpp
	/Users/cjurich/projects/RNAMake/src/structure/beads.cc
	/Users/cjurich/projects/RNAMake/src/structure/residue_type_set.cc
	/Users/cjurich/projects/RNAMake/src/structure/residue_type_set_manager.h 
)
add_library(structure_lib  ${structure_files})
target_link_libraries(structure_lib util_lib  )
###########################################################################
# motif
###########################################################################
set(motif_files
	/Users/cjurich/projects/RNAMake/src/motif/motif_state.cpp
	/Users/cjurich/projects/RNAMake/src/motif/motif_factory.cc
	/Users/cjurich/projects/RNAMake/src/motif/pose_factory.cc
	/Users/cjurich/projects/RNAMake/src/motif/motif_to_secondary_structure.cc
	/Users/cjurich/projects/RNAMake/src/motif/motif_ensemble.cpp
	/Users/cjurich/projects/RNAMake/src/motif/motif_state_ensemble.cpp
	/Users/cjurich/projects/RNAMake/src/motif/motif_state_aligner.cpp
	/Users/cjurich/projects/RNAMake/src/motif/pose.cc
	/Users/cjurich/projects/RNAMake/src/motif/motif_scorer.cc
	/Users/cjurich/projects/RNAMake/src/motif/motif.cc 
)
add_library(motif_lib  ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib  )
###########################################################################
# motif_tools
###########################################################################
set(motif_tools_files
	/Users/cjurich/projects/RNAMake/src/motif_tools/segmenter.cpp 
)
add_library(motif_tools_lib  ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib  )
###########################################################################
# resources
###########################################################################
set(resources_files
	/Users/cjurich/projects/RNAMake/src/resources/resource_manager.cc
	/Users/cjurich/projects/RNAMake/src/resources/motif_state_sqlite_library.cpp
	/Users/cjurich/projects/RNAMake/src/resources/motif_ensemble_sqlite_connection.cpp
	/Users/cjurich/projects/RNAMake/src/resources/motif_state_ensemble_sqlite_library.cpp
	/Users/cjurich/projects/RNAMake/src/resources/added_motif_library.cc
	/Users/cjurich/projects/RNAMake/src/resources/motif_sqlite_connection.cc
	/Users/cjurich/projects/RNAMake/src/resources/sqlite_library.cc
	/Users/cjurich/projects/RNAMake/src/resources/motif_sqlite_library.cc 
)
add_library(resources_lib  ${resources_files})
target_link_libraries(resources_lib motif_lib  )
###########################################################################
# motif_data_structure
###########################################################################
set(motif_data_structure_files
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_graph.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_ensemble_tree.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_node.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_tree.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_topology.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_ensemble_graph.cc
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_connection.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_graph.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_tree.fwd.h
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_merger.cpp
	/Users/cjurich/projects/RNAMake/src/motif_data_structure/motif_state_tree.cpp 
)
add_library(motif_data_structure_lib  ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib  )
###########################################################################
# thermo_fluctuation
###########################################################################
set(thermo_fluctuation_files
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/graph/scorer.cc
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/graph/logging.cc
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation.cpp
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/graph/sampler.cc
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/graph/simulation.cc
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/thermo_fluc_scorer.cpp
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/thermo_fluc_relax.cc
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/thermo_fluc_sampler.cpp
	/Users/cjurich/projects/RNAMake/src/thermo_fluctuation/graph/sterics.cc 
)
add_library(thermo_fluctuation_lib  ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib  )
###########################################################################
# motif_search
###########################################################################
set(motif_search_files
	/Users/cjurich/projects/RNAMake/src/motif_search/monte_carlo/scorer.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/exhaustive/search.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/exhaustive/motif_state_enumerator.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/exhaustive/scorer.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/path_finding/node.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/path_finding/selector.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/monte_carlo/search.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/search.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/solution_filter.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/motif_state_monte_carlo.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/path_finding/scorer.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/path_finding/search.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/problem.cc
	/Users/cjurich/projects/RNAMake/src/motif_search/solution_topology.cc 
)
add_library(motif_search_lib  ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib  )
###########################################################################
# sequence_optimization
###########################################################################
set(sequence_optimization_files
	/Users/cjurich/projects/RNAMake/src/sequence_optimization/sequence_optimizer.cpp
	/Users/cjurich/projects/RNAMake/src/sequence_optimization/sequence_optimizer_3d.cpp 
)
add_library(sequence_optimization_lib  ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib eternabot_lib  )
add_library(all_lib  /Users/cjurich/projects/RNAMake/src/main.cpp)
target_link_libraries(all_lib motif_tools_lib thermo_fluctuation_lib motif_search_lib sequence_optimization_lib )
####################################################################################################
####################################################################################################
# Unittest Declarations 
####################################################################################################
###########################################################################
# base Tests 
###########################################################################
	add_executable(cl_option_unittest /Users/cjurich/projects/RNAMake/unittests/base/cl_option_unittest.cpp)
	target_link_libraries(cl_option_unittest base_lib  )
	add_executable(log_unittest /Users/cjurich/projects/RNAMake/unittests/base/log_unittest.cpp)
	target_link_libraries(log_unittest base_lib  )
	add_executable(string_unittest /Users/cjurich/projects/RNAMake/unittests/base/string_unittest.cpp)
	target_link_libraries(string_unittest base_lib  )
	add_executable(option_unittest /Users/cjurich/projects/RNAMake/unittests/base/option_unittest.cpp)
	target_link_libraries(option_unittest base_lib  )
	add_executable(application_unittest /Users/cjurich/projects/RNAMake/unittests/base/application_unittest.cpp)
	target_link_libraries(application_unittest base_lib  )
	add_executable(env_manager_unittest /Users/cjurich/projects/RNAMake/unittests/base/env_manager_unittest.cpp)
	target_link_libraries(env_manager_unittest base_lib  )
###########################################################################
# math Tests 
###########################################################################
	add_executable(xyz_matrix_unittest /Users/cjurich/projects/RNAMake/unittests/math/xyz_matrix_unittest.cpp)
	target_link_libraries(xyz_matrix_unittest math_lib  )
	add_executable(hashing_unittest /Users/cjurich/projects/RNAMake/unittests/math/hashing_unittest.cpp)
	target_link_libraries(hashing_unittest math_lib  )
	add_executable(quaternion_unittest /Users/cjurich/projects/RNAMake/unittests/math/quaternion_unittest.cpp)
	target_link_libraries(quaternion_unittest math_lib  )
	add_executable(rotation_unittest /Users/cjurich/projects/RNAMake/unittests/math/rotation_unittest.cpp)
	target_link_libraries(rotation_unittest math_lib  )
###########################################################################
# data_structure Tests 
###########################################################################
	add_executable(tree_unittest /Users/cjurich/projects/RNAMake/unittests/data_structure/tree_unittest.cpp)
	target_link_libraries(tree_unittest data_structure_lib  )
	add_executable(graph_iter_list_unittest /Users/cjurich/projects/RNAMake/unittests/data_structure/graph_iter_list_unittest.cpp)
	target_link_libraries(graph_iter_list_unittest data_structure_lib  )
	add_executable(graph_unittest /Users/cjurich/projects/RNAMake/unittests/data_structure/graph_unittest.cpp)
	target_link_libraries(graph_unittest data_structure_lib  )
	add_executable(graph_adjacency_list_unittest /Users/cjurich/projects/RNAMake/unittests/data_structure/graph_adjacency_list_unittest.cpp)
	target_link_libraries(graph_adjacency_list_unittest data_structure_lib  )
###########################################################################
# util Tests 
###########################################################################
	add_executable(steric_lookup_unittest /Users/cjurich/projects/RNAMake/unittests/util/steric_lookup_unittest.cpp)
	target_link_libraries(steric_lookup_unittest util_lib  )
	add_executable(motif_type_unittest /Users/cjurich/projects/RNAMake/unittests/util/motif_type_unittest.cpp)
	target_link_libraries(motif_type_unittest util_lib  )
	add_executable(csv_unittest /Users/cjurich/projects/RNAMake/unittests/util/csv_unittest.cpp)
	target_link_libraries(csv_unittest util_lib  )
	add_executable(x3dna_unittest /Users/cjurich/projects/RNAMake/unittests/util/x3dna_unittest.cpp)
	target_link_libraries(x3dna_unittest util_lib  )
	add_executable(sqlite3_connection_unittest /Users/cjurich/projects/RNAMake/unittests/util/sqlite3_connection_unittest.cpp)
	target_link_libraries(sqlite3_connection_unittest util_lib  )
	add_executable(uuid_unittest /Users/cjurich/projects/RNAMake/unittests/util/uuid_unittest.cpp)
	target_link_libraries(uuid_unittest util_lib  )
	add_executable(dssr_unittest /Users/cjurich/projects/RNAMake/unittests/util/dssr_unittest.cpp)
	target_link_libraries(dssr_unittest util_lib  )
###########################################################################
# vienna Tests 
###########################################################################
###########################################################################
# secondary_structure Tests 
###########################################################################
	add_executable(ss_motif_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/ss_motif_unittest.cpp)
	target_link_libraries(ss_motif_unittest secondary_structure_lib  )
	add_executable(ss_chain_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/ss_chain_unittest.cpp)
	target_link_libraries(ss_chain_unittest secondary_structure_lib  )
	add_executable(ss_residue_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/ss_residue_unittest.cpp)
	target_link_libraries(ss_residue_unittest secondary_structure_lib  )
	add_executable(secondary_structure_parser_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/secondary_structure_parser_unittest.cpp)
	target_link_libraries(secondary_structure_parser_unittest secondary_structure_lib  )
	add_executable(ss_structure_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/ss_structure_unittest.cpp)
	target_link_libraries(ss_structure_unittest secondary_structure_lib  )
	add_executable(ss_pose_unittest /Users/cjurich/projects/RNAMake/unittests/secondary_structure/ss_pose_unittest.cpp)
	target_link_libraries(ss_pose_unittest secondary_structure_lib  )
###########################################################################
# eternabot Tests 
###########################################################################
###########################################################################
# structure Tests 
###########################################################################
	add_executable(atom_unittest /Users/cjurich/projects/RNAMake/unittests/structure/atom_unittest.cpp)
	target_link_libraries(atom_unittest structure_lib  )
	add_executable(residue_unittest /Users/cjurich/projects/RNAMake/unittests/structure/residue_unittest.cpp)
	target_link_libraries(residue_unittest structure_lib  )
	add_executable(basepair_state_unittest /Users/cjurich/projects/RNAMake/unittests/structure/basepair_state_unittest.cpp)
	target_link_libraries(basepair_state_unittest structure_lib  )
	add_executable(basepair_unittest /Users/cjurich/projects/RNAMake/unittests/structure/basepair_unittest.cpp)
	target_link_libraries(basepair_unittest structure_lib  )
	add_executable(chain_closure_unittest /Users/cjurich/projects/RNAMake/unittests/structure/chain_closure_unittest.cpp)
	target_link_libraries(chain_closure_unittest structure_lib  )
	add_executable(residue_type_unittest /Users/cjurich/projects/RNAMake/unittests/structure/residue_type_unittest.cpp)
	target_link_libraries(residue_type_unittest structure_lib  )
	add_executable(chain_unittest /Users/cjurich/projects/RNAMake/unittests/structure/chain_unittest.cpp)
	target_link_libraries(chain_unittest structure_lib  )
	add_executable(pdb_parser_unittest /Users/cjurich/projects/RNAMake/unittests/structure/pdb_parser_unittest.cpp)
	target_link_libraries(pdb_parser_unittest structure_lib  )
	add_executable(structure_unittest /Users/cjurich/projects/RNAMake/unittests/structure/structure_unittest.cpp)
	target_link_libraries(structure_unittest structure_lib  )
	add_executable(cif_parser_unittest /Users/cjurich/projects/RNAMake/unittests/structure/cif_parser_unittest.cpp)
	target_link_libraries(cif_parser_unittest structure_lib  )
###########################################################################
# motif Tests 
###########################################################################
	add_executable(motif_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_unittest.cpp)
	target_link_libraries(motif_unittest motif_lib  )
	add_executable(motif_factory_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_factory_unittest.cpp)
	target_link_libraries(motif_factory_unittest motif_lib  )
	add_executable(motif_state_ensemble_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_state_ensemble_unittest.cpp)
	target_link_libraries(motif_state_ensemble_unittest motif_lib  )
	add_executable(motif_scorer_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_scorer_unittest.cpp)
	target_link_libraries(motif_scorer_unittest motif_lib  )
	add_executable(motif_state_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_state_unittest.cpp)
	target_link_libraries(motif_state_unittest motif_lib  )
	add_executable(motif_to_secondary_structure_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_to_secondary_structure_unittest.cpp)
	target_link_libraries(motif_to_secondary_structure_unittest motif_lib  )
	add_executable(motif_ensemble_unittest /Users/cjurich/projects/RNAMake/unittests/motif/motif_ensemble_unittest.cpp)
	target_link_libraries(motif_ensemble_unittest motif_lib  )
###########################################################################
# motif_tools Tests 
###########################################################################
	add_executable(segmenter_unittest /Users/cjurich/projects/RNAMake/unittests/motif_tools/segmenter_unittest.cpp)
	target_link_libraries(segmenter_unittest motif_tools_lib  )
###########################################################################
# resources Tests 
###########################################################################
	add_executable(added_motif_library_unittest /Users/cjurich/projects/RNAMake/unittests/resources/added_motif_library_unittest.cpp)
	target_link_libraries(added_motif_library_unittest resources_lib  )
	add_executable(motif_state_ensemble_sqlite_library_unittest /Users/cjurich/projects/RNAMake/unittests/resources/motif_state_ensemble_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_ensemble_sqlite_library_unittest resources_lib  )
	add_executable(motif_ensemble_sqlite_connection_unittest /Users/cjurich/projects/RNAMake/unittests/resources/motif_ensemble_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_ensemble_sqlite_connection_unittest resources_lib  )
	add_executable(resource_mananger_unittest /Users/cjurich/projects/RNAMake/unittests/resources/resource_mananger_unittest.cpp)
	target_link_libraries(resource_mananger_unittest resources_lib  )
	add_executable(motif_state_sqlite_library_unittest /Users/cjurich/projects/RNAMake/unittests/resources/motif_state_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_sqlite_library_unittest resources_lib  )
	add_executable(motif_sqlite_library_unittest /Users/cjurich/projects/RNAMake/unittests/resources/motif_sqlite_library_unittest.cpp)
	target_link_libraries(motif_sqlite_library_unittest resources_lib  )
	add_executable(motif_sqlite_connection_unittest /Users/cjurich/projects/RNAMake/unittests/resources/motif_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_sqlite_connection_unittest resources_lib  )
###########################################################################
# motif_data_structure Tests 
###########################################################################
	add_executable(motif_graph_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_graph_unittest.cpp)
	target_link_libraries(motif_graph_unittest motif_data_structure_lib  )
	add_executable(motif_state_ensemble_tree_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_state_ensemble_tree_unittest.cpp)
	target_link_libraries(motif_state_ensemble_tree_unittest motif_data_structure_lib  )
	add_executable(motif_merger_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_merger_unittest.cpp)
	target_link_libraries(motif_merger_unittest motif_data_structure_lib  )
	add_executable(motif_tree_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_tree_unittest.cpp)
	target_link_libraries(motif_tree_unittest motif_data_structure_lib  )
	add_executable(motif_state_tree_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_state_tree_unittest.cpp)
	target_link_libraries(motif_state_tree_unittest motif_data_structure_lib  )
	add_executable(motif_topology_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_topology_unittest.cpp)
	target_link_libraries(motif_topology_unittest motif_data_structure_lib  )
	add_executable(motif_state_graph_unittest /Users/cjurich/projects/RNAMake/unittests/motif_data_structure/motif_state_graph_unittest.cpp)
	target_link_libraries(motif_state_graph_unittest motif_data_structure_lib  )
###########################################################################
# thermo_fluctuation Tests 
###########################################################################
	add_executable(thermo_fluc_simulation_unittest /Users/cjurich/projects/RNAMake/unittests/thermo_fluctuation/thermo_fluc_simulation_unittest.cpp)
	target_link_libraries(thermo_fluc_simulation_unittest thermo_fluctuation_lib  )
	add_executable(thermo_fluc_sampler_unittest /Users/cjurich/projects/RNAMake/unittests/thermo_fluctuation/thermo_fluc_sampler_unittest.cpp)
	target_link_libraries(thermo_fluc_sampler_unittest thermo_fluctuation_lib  )
###########################################################################
# motif_search Tests 
###########################################################################
	add_executable(monte_carlo_search_unittest /Users/cjurich/projects/RNAMake/unittests/motif_search/monte_carlo_search_unittest.cpp)
	target_link_libraries(monte_carlo_search_unittest motif_search_lib  )
	add_executable(exhaustive_search_unittest /Users/cjurich/projects/RNAMake/unittests/motif_search/exhaustive_search_unittest.cpp)
	target_link_libraries(exhaustive_search_unittest motif_search_lib  )
	add_executable(motif_state_search_unittest /Users/cjurich/projects/RNAMake/unittests/motif_search/motif_state_search_unittest.cpp)
	target_link_libraries(motif_state_search_unittest motif_search_lib  )
###########################################################################
# sequence_optimization Tests 
###########################################################################
	add_executable(sequence_optimizer_unittest /Users/cjurich/projects/RNAMake/unittests/sequence_optimization/sequence_optimizer_unittest.cpp)
	target_link_libraries(sequence_optimizer_unittest sequence_optimization_lib  )
####################################################################################################
####################################################################################################
# App declarations
###############################################################################################################################################################################
# ddg_tecto
###########################################################################
	add_executable(ddg_tecto /Users/cjurich/projects/RNAMake/apps/ddg_tecto/ddg_tecto.cc)
	target_link_libraries(ddg_tecto all_lib   )
###########################################################################
# build_sqlite_libraries
###########################################################################
	add_executable(build_sqlite_libraries /Users/cjurich/projects/RNAMake/apps/build_sqlite_libraries/build_sqlite_libraries.cpp)
	target_link_libraries(build_sqlite_libraries all_lib   )
###########################################################################
# design_rna_scaffold
###########################################################################
	add_executable(design_rna_scaffold /Users/cjurich/projects/RNAMake/apps/design_rna_scaffold/design_rna_scaffold.cc)
	target_link_libraries(design_rna_scaffold all_lib   )
####################################################################################################
