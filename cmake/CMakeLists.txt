####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
project(RNAMake)

set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed ")
set(CMAKE_EXE_LINKER_FLAGS " -lstdc++ -Wl,--no-as-needed,--no-export-dynamic ")
include(/home/hassan/RNAMake/cmake/build/compiler.cmake)
set(RNAMAKE /home/hassan/RNAMake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
include_directories(/home/hassan/RNAMake/src/)
include_directories(/home/hassan/RNAMake/src/external/)
include_directories(/home/hassan/RNAMake/unittests/)
include_directories(/home/hassan/RNAMake/apps/)
include(/home/hassan/RNAMake/cmake/build/sqlite.cmake)
 ####################################################################################################
####################################################################################################
# Library declarations 
####################################################################################################
###########################################################################
# base
###########################################################################
set(base_files
	/home/hassan/RNAMake/src/base/log.cc
	/home/hassan/RNAMake/src/base/types.h
	/home/hassan/RNAMake/src/base/option.cc
	/home/hassan/RNAMake/src/base/util.hpp
	/home/hassan/RNAMake/src/base/string.cc
	/home/hassan/RNAMake/src/base/vector_container.h
	/home/hassan/RNAMake/src/base/file_io.cc
	/home/hassan/RNAMake/src/base/application.cpp
	/home/hassan/RNAMake/src/base/settings.cc
	/home/hassan/RNAMake/src/base/backtrace.cpp
	/home/hassan/RNAMake/src/base/command_line_parser.cpp
	/home/hassan/RNAMake/src/base/global_constants.cc
	/home/hassan/RNAMake/src/base/env_manager.cpp
	/home/hassan/RNAMake/src/base/assertions.h
	/home/hassan/RNAMake/src/base/exception.cpp
	/home/hassan/RNAMake/src/base/cl_option.cc
	/home/hassan/RNAMake/src/base/sys_interface.cpp 
)
add_library(base_lib STATIC ${base_files})
target_link_libraries(base_lib  -static)
###########################################################################
# math
###########################################################################
set(math_files
	/home/hassan/RNAMake/src/math/numerical.cc
	/home/hassan/RNAMake/src/math/xyz_vector.h
	/home/hassan/RNAMake/src/math/euler.cc
	/home/hassan/RNAMake/src/math/transform.h
	/home/hassan/RNAMake/src/math/xyz_matrix.h
	/home/hassan/RNAMake/src/math/stats.cc
	/home/hassan/RNAMake/src/math/hashing.cc
	/home/hassan/RNAMake/src/math/quaternion.cc 
)
add_library(math_lib STATIC ${math_files})
target_link_libraries(math_lib base_lib  -static)
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	/home/hassan/RNAMake/src/data_structure/tree/tree_node.cc
	/home/hassan/RNAMake/src/data_structure/graph/graph.cc
	/home/hassan/RNAMake/src/data_structure/tree/tree.cc
	/home/hassan/RNAMake/src/data_structure/graph_base.cc
	/home/hassan/RNAMake/src/data_structure/graph_iter_list.cc
	/home/hassan/RNAMake/src/data_structure/graph.cc
	/home/hassan/RNAMake/src/data_structure/graph/graph_node.fwd.h
	/home/hassan/RNAMake/src/data_structure/graph/graph_node.cc
	/home/hassan/RNAMake/src/data_structure/graph_adjacency_list.cc 
)
add_library(data_structure_lib STATIC ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  -static)
###########################################################################
# util
###########################################################################
set(util_files
	/home/hassan/RNAMake/src/util/cartesian_product.cc
	/home/hassan/RNAMake/src/util/basic_io.cpp
	/home/hassan/RNAMake/src/util/sqlite3_connection.cc
	/home/hassan/RNAMake/src/util/random_number_generator.h
	/home/hassan/RNAMake/src/util/motif_type.cc
	/home/hassan/RNAMake/src/util/csv.h
	/home/hassan/RNAMake/src/util/uuid.cc
	/home/hassan/RNAMake/src/util/dssr.cpp
	/home/hassan/RNAMake/src/util/monte_carlo.cpp
	/home/hassan/RNAMake/src/util/steric_lookup.cpp
	/home/hassan/RNAMake/src/util/x3dna.cc 
)
add_library(util_lib STATIC ${util_files})
target_link_libraries(util_lib math_lib sqlite3 -static)
###########################################################################
# vienna
###########################################################################
set(vienna_files
	/home/hassan/RNAMake/src/vienna/energy_par.cc
	/home/hassan/RNAMake/src/vienna/vienna.cc
	/home/hassan/RNAMake/src/vienna/pair_mat.h 
)
add_library(vienna_lib STATIC ${vienna_files})
target_link_libraries(vienna_lib base_lib  -static)
###########################################################################
# secondary_structure
###########################################################################
set(secondary_structure_files
	/home/hassan/RNAMake/src/secondary_structure/sequence_constraint.cc
	/home/hassan/RNAMake/src/secondary_structure/util.cpp
	/home/hassan/RNAMake/src/secondary_structure/structure.cpp
	/home/hassan/RNAMake/src/secondary_structure/pose.cpp
	/home/hassan/RNAMake/src/secondary_structure/secondary_structure_tree.cpp
	/home/hassan/RNAMake/src/secondary_structure/basepair.cc
	/home/hassan/RNAMake/src/secondary_structure/sequence_tools.cc
	/home/hassan/RNAMake/src/secondary_structure/chain.cc
	/home/hassan/RNAMake/src/secondary_structure/secondary_structure_parser.cpp
	/home/hassan/RNAMake/src/secondary_structure/motif.cc
	/home/hassan/RNAMake/src/secondary_structure/rna_structure.cpp
	/home/hassan/RNAMake/src/secondary_structure/residue.cc 
)
add_library(secondary_structure_lib STATIC ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib  -static)
###########################################################################
# eternabot
###########################################################################
set(eternabot_files
	/home/hassan/RNAMake/src/eternabot/strategy/a_basic_test.cpp
	/home/hassan/RNAMake/src/eternabot/strategy/berex_test.cpp
	/home/hassan/RNAMake/src/eternabot/strategy/clear_plot.cpp
	/home/hassan/RNAMake/src/eternabot/feature_generator.cpp
	/home/hassan/RNAMake/src/eternabot/strategy/num_of_yellow.cpp
	/home/hassan/RNAMake/src/eternabot/scorer.cpp
	/home/hassan/RNAMake/src/eternabot/strategy.cpp
	/home/hassan/RNAMake/src/eternabot/strategy/direction_of_gc.cpp
	/home/hassan/RNAMake/src/eternabot/sequence_designer.cpp 
)
add_library(eternabot_lib STATIC ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib secondary_structure_lib  -static)
###########################################################################
# structure
###########################################################################
set(structure_files
	/home/hassan/RNAMake/src/structure/is_equal.cpp
	/home/hassan/RNAMake/src/structure/residue_type_set_manager.h
	/home/hassan/RNAMake/src/structure/structure.cc
	/home/hassan/RNAMake/src/structure/residue_type_set.cc
	/home/hassan/RNAMake/src/structure/close_chain.cc
	/home/hassan/RNAMake/src/structure/rna_structure.cpp
	/home/hassan/RNAMake/src/structure/residue.cc
	/home/hassan/RNAMake/src/structure/chain.cc
	/home/hassan/RNAMake/src/structure/residue_type.cc
	/home/hassan/RNAMake/src/structure/cif_parser.cc
	/home/hassan/RNAMake/src/structure/basepair_state.cc
	/home/hassan/RNAMake/src/structure/chain.fwd.h
	/home/hassan/RNAMake/src/structure/pdb_parser.cc
	/home/hassan/RNAMake/src/structure/atom.cc
	/home/hassan/RNAMake/src/structure/beads.cc
	/home/hassan/RNAMake/src/structure/exceptions.cc
	/home/hassan/RNAMake/src/structure/basepair.cc 
)
add_library(structure_lib STATIC ${structure_files})
target_link_libraries(structure_lib util_lib  -static)
###########################################################################
# motif
###########################################################################
set(motif_files
	/home/hassan/RNAMake/src/motif/motif_scorer.cc
	/home/hassan/RNAMake/src/motif/motif_ensemble.cpp
	/home/hassan/RNAMake/src/motif/motif_state.cpp
	/home/hassan/RNAMake/src/motif/motif_state_aligner.cpp
	/home/hassan/RNAMake/src/motif/motif_state_ensemble.cpp
	/home/hassan/RNAMake/src/motif/motif.cc
	/home/hassan/RNAMake/src/motif/motif_factory.cc
	/home/hassan/RNAMake/src/motif/motif_to_secondary_structure.cc
	/home/hassan/RNAMake/src/motif/pose_factory.cc
	/home/hassan/RNAMake/src/motif/pose.cc 
)
add_library(motif_lib STATIC ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib  -static)
###########################################################################
# motif_tools
###########################################################################
set(motif_tools_files
	/home/hassan/RNAMake/src/motif_tools/segmenter.cpp 
)
add_library(motif_tools_lib STATIC ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib  -static)
###########################################################################
# resources
###########################################################################
set(resources_files
	/home/hassan/RNAMake/src/resources/motif_sqlite_library.cc
	/home/hassan/RNAMake/src/resources/motif_ensemble_sqlite_connection.cpp
	/home/hassan/RNAMake/src/resources/motif_sqlite_connection.cc
	/home/hassan/RNAMake/src/resources/resource_manager.cc
	/home/hassan/RNAMake/src/resources/sqlite_library.cc
	/home/hassan/RNAMake/src/resources/motif_state_ensemble_sqlite_library.cpp
	/home/hassan/RNAMake/src/resources/motif_state_sqlite_library.cpp
	/home/hassan/RNAMake/src/resources/added_motif_library.cc 
)
add_library(resources_lib STATIC ${resources_files})
target_link_libraries(resources_lib motif_lib  -static)
###########################################################################
# motif_data_structure
###########################################################################
set(motif_data_structure_files
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_tree.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_tree.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_connection.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_graph.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_ensemble_graph.cc
	/home/hassan/RNAMake/src/motif_data_structure/motif_graph.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_tree.fwd.h
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_ensemble_tree.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_state_node.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_merger.cpp
	/home/hassan/RNAMake/src/motif_data_structure/motif_topology.cpp 
)
add_library(motif_data_structure_lib STATIC ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib  -static)
###########################################################################
# thermo_fluctuation
###########################################################################
set(thermo_fluctuation_files
	/home/hassan/RNAMake/src/thermo_fluctuation/thermo_fluc_sampler.cpp
	/home/hassan/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	/home/hassan/RNAMake/src/thermo_fluctuation/graph/simulation.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/graph/logging.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/thermo_fluc_relax.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/graph/sterics.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/graph/sampler.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/graph/scorer.cc
	/home/hassan/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation.cpp
	/home/hassan/RNAMake/src/thermo_fluctuation/thermo_fluc_scorer.cpp 
)
add_library(thermo_fluctuation_lib STATIC ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib  -static)
###########################################################################
# motif_search
###########################################################################
set(motif_search_files
	/home/hassan/RNAMake/src/motif_search/monte_carlo/search.cc
	/home/hassan/RNAMake/src/motif_search/solution_topology.cc
	/home/hassan/RNAMake/src/motif_search/monte_carlo/scorer.cc
	/home/hassan/RNAMake/src/motif_search/path_finding/selector.cc
	/home/hassan/RNAMake/src/motif_search/problem.cc
	/home/hassan/RNAMake/src/motif_search/exhaustive/motif_state_enumerator.cc
	/home/hassan/RNAMake/src/motif_search/exhaustive/search.cc
	/home/hassan/RNAMake/src/motif_search/path_finding/search.cc
	/home/hassan/RNAMake/src/motif_search/search.cc
	/home/hassan/RNAMake/src/motif_search/exhaustive/scorer.cc
	/home/hassan/RNAMake/src/motif_search/path_finding/scorer.cc
	/home/hassan/RNAMake/src/motif_search/solution_filter.cc
	/home/hassan/RNAMake/src/motif_search/motif_state_monte_carlo.cc
	/home/hassan/RNAMake/src/motif_search/path_finding/node.cc 
)
add_library(motif_search_lib STATIC ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib  -static)
###########################################################################
# sequence_optimization
###########################################################################
set(sequence_optimization_files
	/home/hassan/RNAMake/src/sequence_optimization/sequence_optimizer_3d.cpp
	/home/hassan/RNAMake/src/sequence_optimization/sequence_optimizer.cpp 
)
add_library(sequence_optimization_lib STATIC ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib eternabot_lib  -static)
add_library(all_lib STATIC /home/hassan/RNAMake/src/main.cpp)
target_link_libraries(all_lib motif_tools_lib thermo_fluctuation_lib motif_search_lib sequence_optimization_lib -static)
####################################################################################################
####################################################################################################
# Unittest Declarations 
####################################################################################################
###########################################################################
# base Tests 
###########################################################################
	add_executable(env_manager_unittest /home/hassan/RNAMake/unittests/base/env_manager_unittest.cpp)
	target_link_libraries(env_manager_unittest base_lib -static )
	add_executable(application_unittest /home/hassan/RNAMake/unittests/base/application_unittest.cpp)
	target_link_libraries(application_unittest base_lib -static )
	add_executable(option_unittest /home/hassan/RNAMake/unittests/base/option_unittest.cpp)
	target_link_libraries(option_unittest base_lib -static )
	add_executable(cl_option_unittest /home/hassan/RNAMake/unittests/base/cl_option_unittest.cpp)
	target_link_libraries(cl_option_unittest base_lib -static )
	add_executable(string_unittest /home/hassan/RNAMake/unittests/base/string_unittest.cpp)
	target_link_libraries(string_unittest base_lib -static )
	add_executable(log_unittest /home/hassan/RNAMake/unittests/base/log_unittest.cpp)
	target_link_libraries(log_unittest base_lib -static )
###########################################################################
# math Tests 
###########################################################################
	add_executable(rotation_unittest /home/hassan/RNAMake/unittests/math/rotation_unittest.cpp)
	target_link_libraries(rotation_unittest math_lib -static )
	add_executable(quaternion_unittest /home/hassan/RNAMake/unittests/math/quaternion_unittest.cpp)
	target_link_libraries(quaternion_unittest math_lib -static )
	add_executable(hashing_unittest /home/hassan/RNAMake/unittests/math/hashing_unittest.cpp)
	target_link_libraries(hashing_unittest math_lib -static )
	add_executable(xyz_matrix_unittest /home/hassan/RNAMake/unittests/math/xyz_matrix_unittest.cpp)
	target_link_libraries(xyz_matrix_unittest math_lib -static )
###########################################################################
# data_structure Tests 
###########################################################################
	add_executable(graph_adjacency_list_unittest /home/hassan/RNAMake/unittests/data_structure/graph_adjacency_list_unittest.cpp)
	target_link_libraries(graph_adjacency_list_unittest data_structure_lib -static )
	add_executable(graph_unittest /home/hassan/RNAMake/unittests/data_structure/graph_unittest.cpp)
	target_link_libraries(graph_unittest data_structure_lib -static )
	add_executable(tree_unittest /home/hassan/RNAMake/unittests/data_structure/tree_unittest.cpp)
	target_link_libraries(tree_unittest data_structure_lib -static )
	add_executable(graph_iter_list_unittest /home/hassan/RNAMake/unittests/data_structure/graph_iter_list_unittest.cpp)
	target_link_libraries(graph_iter_list_unittest data_structure_lib -static )
###########################################################################
# util Tests 
###########################################################################
	add_executable(motif_type_unittest /home/hassan/RNAMake/unittests/util/motif_type_unittest.cpp)
	target_link_libraries(motif_type_unittest util_lib -static )
	add_executable(uuid_unittest /home/hassan/RNAMake/unittests/util/uuid_unittest.cpp)
	target_link_libraries(uuid_unittest util_lib -static )
	add_executable(dssr_unittest /home/hassan/RNAMake/unittests/util/dssr_unittest.cpp)
	target_link_libraries(dssr_unittest util_lib -static )
	add_executable(x3dna_unittest /home/hassan/RNAMake/unittests/util/x3dna_unittest.cpp)
	target_link_libraries(x3dna_unittest util_lib -static )
	add_executable(sqlite3_connection_unittest /home/hassan/RNAMake/unittests/util/sqlite3_connection_unittest.cpp)
	target_link_libraries(sqlite3_connection_unittest util_lib -static )
	add_executable(steric_lookup_unittest /home/hassan/RNAMake/unittests/util/steric_lookup_unittest.cpp)
	target_link_libraries(steric_lookup_unittest util_lib -static )
	add_executable(csv_unittest /home/hassan/RNAMake/unittests/util/csv_unittest.cpp)
	target_link_libraries(csv_unittest util_lib -static )
###########################################################################
# vienna Tests 
###########################################################################
###########################################################################
# secondary_structure Tests 
###########################################################################
	add_executable(ss_chain_unittest /home/hassan/RNAMake/unittests/secondary_structure/ss_chain_unittest.cpp)
	target_link_libraries(ss_chain_unittest secondary_structure_lib -static )
	add_executable(ss_structure_unittest /home/hassan/RNAMake/unittests/secondary_structure/ss_structure_unittest.cpp)
	target_link_libraries(ss_structure_unittest secondary_structure_lib -static )
	add_executable(ss_residue_unittest /home/hassan/RNAMake/unittests/secondary_structure/ss_residue_unittest.cpp)
	target_link_libraries(ss_residue_unittest secondary_structure_lib -static )
	add_executable(secondary_structure_parser_unittest /home/hassan/RNAMake/unittests/secondary_structure/secondary_structure_parser_unittest.cpp)
	target_link_libraries(secondary_structure_parser_unittest secondary_structure_lib -static )
	add_executable(ss_pose_unittest /home/hassan/RNAMake/unittests/secondary_structure/ss_pose_unittest.cpp)
	target_link_libraries(ss_pose_unittest secondary_structure_lib -static )
	add_executable(ss_motif_unittest /home/hassan/RNAMake/unittests/secondary_structure/ss_motif_unittest.cpp)
	target_link_libraries(ss_motif_unittest secondary_structure_lib -static )
###########################################################################
# eternabot Tests 
###########################################################################
###########################################################################
# structure Tests 
###########################################################################
	add_executable(basepair_state_unittest /home/hassan/RNAMake/unittests/structure/basepair_state_unittest.cpp)
	target_link_libraries(basepair_state_unittest structure_lib -static )
	add_executable(residue_unittest /home/hassan/RNAMake/unittests/structure/residue_unittest.cpp)
	target_link_libraries(residue_unittest structure_lib -static )
	add_executable(basepair_unittest /home/hassan/RNAMake/unittests/structure/basepair_unittest.cpp)
	target_link_libraries(basepair_unittest structure_lib -static )
	add_executable(pdb_parser_unittest /home/hassan/RNAMake/unittests/structure/pdb_parser_unittest.cpp)
	target_link_libraries(pdb_parser_unittest structure_lib -static )
	add_executable(cif_parser_unittest /home/hassan/RNAMake/unittests/structure/cif_parser_unittest.cpp)
	target_link_libraries(cif_parser_unittest structure_lib -static )
	add_executable(atom_unittest /home/hassan/RNAMake/unittests/structure/atom_unittest.cpp)
	target_link_libraries(atom_unittest structure_lib -static )
	add_executable(chain_unittest /home/hassan/RNAMake/unittests/structure/chain_unittest.cpp)
	target_link_libraries(chain_unittest structure_lib -static )
	add_executable(structure_unittest /home/hassan/RNAMake/unittests/structure/structure_unittest.cpp)
	target_link_libraries(structure_unittest structure_lib -static )
	add_executable(chain_closure_unittest /home/hassan/RNAMake/unittests/structure/chain_closure_unittest.cpp)
	target_link_libraries(chain_closure_unittest structure_lib -static )
	add_executable(residue_type_unittest /home/hassan/RNAMake/unittests/structure/residue_type_unittest.cpp)
	target_link_libraries(residue_type_unittest structure_lib -static )
###########################################################################
# motif Tests 
###########################################################################
	add_executable(motif_ensemble_unittest /home/hassan/RNAMake/unittests/motif/motif_ensemble_unittest.cpp)
	target_link_libraries(motif_ensemble_unittest motif_lib -static )
	add_executable(motif_unittest /home/hassan/RNAMake/unittests/motif/motif_unittest.cpp)
	target_link_libraries(motif_unittest motif_lib -static )
	add_executable(motif_state_ensemble_unittest /home/hassan/RNAMake/unittests/motif/motif_state_ensemble_unittest.cpp)
	target_link_libraries(motif_state_ensemble_unittest motif_lib -static )
	add_executable(motif_factory_unittest /home/hassan/RNAMake/unittests/motif/motif_factory_unittest.cpp)
	target_link_libraries(motif_factory_unittest motif_lib -static )
	add_executable(motif_to_secondary_structure_unittest /home/hassan/RNAMake/unittests/motif/motif_to_secondary_structure_unittest.cpp)
	target_link_libraries(motif_to_secondary_structure_unittest motif_lib -static )
	add_executable(motif_scorer_unittest /home/hassan/RNAMake/unittests/motif/motif_scorer_unittest.cpp)
	target_link_libraries(motif_scorer_unittest motif_lib -static )
	add_executable(motif_state_unittest /home/hassan/RNAMake/unittests/motif/motif_state_unittest.cpp)
	target_link_libraries(motif_state_unittest motif_lib -static )
###########################################################################
# motif_tools Tests 
###########################################################################
	add_executable(segmenter_unittest /home/hassan/RNAMake/unittests/motif_tools/segmenter_unittest.cpp)
	target_link_libraries(segmenter_unittest motif_tools_lib -static )
###########################################################################
# resources Tests 
###########################################################################
	add_executable(motif_state_ensemble_sqlite_library_unittest /home/hassan/RNAMake/unittests/resources/motif_state_ensemble_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_ensemble_sqlite_library_unittest resources_lib -static )
	add_executable(motif_sqlite_library_unittest /home/hassan/RNAMake/unittests/resources/motif_sqlite_library_unittest.cpp)
	target_link_libraries(motif_sqlite_library_unittest resources_lib -static )
	add_executable(motif_sqlite_connection_unittest /home/hassan/RNAMake/unittests/resources/motif_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_sqlite_connection_unittest resources_lib -static )
	add_executable(resource_mananger_unittest /home/hassan/RNAMake/unittests/resources/resource_mananger_unittest.cpp)
	target_link_libraries(resource_mananger_unittest resources_lib -static )
	add_executable(motif_state_sqlite_library_unittest /home/hassan/RNAMake/unittests/resources/motif_state_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_sqlite_library_unittest resources_lib -static )
	add_executable(added_motif_library_unittest /home/hassan/RNAMake/unittests/resources/added_motif_library_unittest.cpp)
	target_link_libraries(added_motif_library_unittest resources_lib -static )
	add_executable(motif_ensemble_sqlite_connection_unittest /home/hassan/RNAMake/unittests/resources/motif_ensemble_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_ensemble_sqlite_connection_unittest resources_lib -static )
###########################################################################
# motif_data_structure Tests 
###########################################################################
	add_executable(motif_graph_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_graph_unittest.cpp)
	target_link_libraries(motif_graph_unittest motif_data_structure_lib -static )
	add_executable(motif_state_graph_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_state_graph_unittest.cpp)
	target_link_libraries(motif_state_graph_unittest motif_data_structure_lib -static )
	add_executable(motif_state_tree_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_state_tree_unittest.cpp)
	target_link_libraries(motif_state_tree_unittest motif_data_structure_lib -static )
	add_executable(motif_topology_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_topology_unittest.cpp)
	target_link_libraries(motif_topology_unittest motif_data_structure_lib -static )
	add_executable(motif_merger_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_merger_unittest.cpp)
	target_link_libraries(motif_merger_unittest motif_data_structure_lib -static )
	add_executable(motif_tree_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_tree_unittest.cpp)
	target_link_libraries(motif_tree_unittest motif_data_structure_lib -static )
	add_executable(motif_state_ensemble_tree_unittest /home/hassan/RNAMake/unittests/motif_data_structure/motif_state_ensemble_tree_unittest.cpp)
	target_link_libraries(motif_state_ensemble_tree_unittest motif_data_structure_lib -static )
###########################################################################
# thermo_fluctuation Tests 
###########################################################################
	add_executable(thermo_fluc_simulation_unittest /home/hassan/RNAMake/unittests/thermo_fluctuation/thermo_fluc_simulation_unittest.cpp)
	target_link_libraries(thermo_fluc_simulation_unittest thermo_fluctuation_lib -static )
	add_executable(thermo_fluc_sampler_unittest /home/hassan/RNAMake/unittests/thermo_fluctuation/thermo_fluc_sampler_unittest.cpp)
	target_link_libraries(thermo_fluc_sampler_unittest thermo_fluctuation_lib -static )
###########################################################################
# motif_search Tests 
###########################################################################
	add_executable(motif_state_search_unittest /home/hassan/RNAMake/unittests/motif_search/motif_state_search_unittest.cpp)
	target_link_libraries(motif_state_search_unittest motif_search_lib -static )
	add_executable(monte_carlo_search_unittest /home/hassan/RNAMake/unittests/motif_search/monte_carlo_search_unittest.cpp)
	target_link_libraries(monte_carlo_search_unittest motif_search_lib -static )
	add_executable(exhaustive_search_unittest /home/hassan/RNAMake/unittests/motif_search/exhaustive_search_unittest.cpp)
	target_link_libraries(exhaustive_search_unittest motif_search_lib -static )
###########################################################################
# sequence_optimization Tests 
###########################################################################
	add_executable(sequence_optimizer_unittest /home/hassan/RNAMake/unittests/sequence_optimization/sequence_optimizer_unittest.cpp)
	target_link_libraries(sequence_optimizer_unittest sequence_optimization_lib -static )
####################################################################################################
####################################################################################################
# App declarations
####################################################################################################
###########################################################################
# ddg_tecto
###########################################################################
	add_executable( ddg_tecto /home/hassan/RNAMake/apps/ddg_tecto/ddg_tecto.cc)
	target_link_libraries(ddg_tecto all_lib -static  )
###########################################################################
# tecto_graph_setup
###########################################################################
	add_executable( tecto_graph_setup /home/hassan/RNAMake/apps/ddg_tecto/tecto_graph_setup.cc)
	target_link_libraries(tecto_graph_setup all_lib -static  )
###########################################################################
# build_sqlite_libraries
###########################################################################
	add_executable( build_sqlite_libraries /home/hassan/RNAMake/apps/build_sqlite_libraries/build_sqlite_libraries.cpp)
	target_link_libraries(build_sqlite_libraries all_lib -static  )
###########################################################################
# design_rna_scaffold
###########################################################################
	add_executable( design_rna_scaffold /home/hassan/RNAMake/apps/design_rna_scaffold/design_rna_scaffold.cc)
	target_link_libraries(design_rna_scaffold all_lib -static  )
###########################################################################
# build_motif_graph
###########################################################################
	add_executable( build_motif_graph /home/hassan/RNAMake/apps/build_motif_graph/build_motif_graph.cxx)
	target_link_libraries(build_motif_graph all_lib -static  )
###########################################################################
# compare_db_files
###########################################################################
	add_executable( compare_db_files /home/hassan/RNAMake/apps/compare_db_files/compare_db_files.cpp)
	target_link_libraries(compare_db_files all_lib -static  )
####################################################################################################
