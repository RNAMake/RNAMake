####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
project(rnamake_new)

set(CMAKE_CXX_STANDARD 14)
set( CMAKE_CXX_FLAGS " -pthread -L/opt/local/lib " )
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-as-needed -ldl")
set(CMAKE_EXE_LINKER_FLAGS " -lstdc++ -Wl,--no-as-needed -ldl ")
include(/RNAMake/cmake/build/compiler.cmake)

# Include path for Python header files
# Include paths for RNAMake src
include_directories(/RNAMake/src/)
include_directories(/RNAMake/src/plog/)
include_directories(/RNAMake/unittests/)

include_directories(/RNAMake/apps/)

# sqlite libraries
find_library(SQLITE3_LIBRARY sqlite3 VARIANT static )
####################################################################################################
####################################################################################################
# Library declarations 
####################################################################################################
###########################################################################
# base
###########################################################################
set(base_files
	/RNAMake/src/base/file_io.cc
	/RNAMake/src/base/util.hpp
	/RNAMake/src/base/assertions.h
	/RNAMake/src/base/backtrace.cpp
	/RNAMake/src/base/log.cc
	/RNAMake/src/base/types.h
	/RNAMake/src/base/cl_option.cc
	/RNAMake/src/base/string.cc
	/RNAMake/src/base/option.cc
	/RNAMake/src/base/vector_container.h
	/RNAMake/src/base/global_constants.cc
	/RNAMake/src/base/settings.cc
	/RNAMake/src/base/application.cpp
	/RNAMake/src/base/command_line_parser.cpp 
)
add_library(base_lib STATIC ${base_files})
target_link_libraries(base_lib  -static)
###########################################################################
# math
###########################################################################
set(math_files
	/RNAMake/src/math/euler.cc
	/RNAMake/src/math/xyz_matrix.h
	/RNAMake/src/math/quaternion.cc
	/RNAMake/src/math/transform.h
	/RNAMake/src/math/hashing.cc
	/RNAMake/src/math/numerical.cc
	/RNAMake/src/math/stats.cc
	/RNAMake/src/math/xyz_vector.h 
)
add_library(math_lib STATIC ${math_files})
target_link_libraries(math_lib base_lib  -static)
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	/RNAMake/src/data_structure/graph_iter_list.cc
	/RNAMake/src/data_structure/tree/tree.cc
	/RNAMake/src/data_structure/graph.cc
	/RNAMake/src/data_structure/tree/tree_node.cc
	/RNAMake/src/data_structure/graph_base.cc
	/RNAMake/src/data_structure/graph/graph.cc
	/RNAMake/src/data_structure/graph/graph_node.cc
	/RNAMake/src/data_structure/graph/graph_node.fwd.h
	/RNAMake/src/data_structure/graph_adjacency_list.cc 
)
add_library(data_structure_lib STATIC ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  -static)
###########################################################################
# util
###########################################################################
set(util_files
	/RNAMake/src/util/monte_carlo.cpp
	/RNAMake/src/util/motif_type.cc
	/RNAMake/src/util/sqlite3_connection.cc
	/RNAMake/src/util/csv.cpp
	/RNAMake/src/util/x3dna.cc
	/RNAMake/src/util/uuid.cc
	/RNAMake/src/util/cartesian_product.cc
	/RNAMake/src/util/random_number_generator.h
	/RNAMake/src/util/basic_io.cpp
	/RNAMake/src/util/steric_lookup.cpp 
)
add_library(util_lib STATIC ${util_files})
target_link_libraries(util_lib math_lib ${SQLITE3_LIBRARY} "-ldl" -static)
###########################################################################
# vienna
###########################################################################
set(vienna_files
	/RNAMake/src/vienna/pair_mat.h
	/RNAMake/src/vienna/vienna.cc
	/RNAMake/src/vienna/energy_par.cc 
)
add_library(vienna_lib STATIC ${vienna_files})
target_link_libraries(vienna_lib base_lib  -static)
###########################################################################
# secondary_structure
###########################################################################
set(secondary_structure_files
	/RNAMake/src/secondary_structure/sequence_tools.cc
	/RNAMake/src/secondary_structure/util.cpp
	/RNAMake/src/secondary_structure/secondary_structure_parser.cpp
	/RNAMake/src/secondary_structure/chain.cc
	/RNAMake/src/secondary_structure/motif.cc
	/RNAMake/src/secondary_structure/structure.cpp
	/RNAMake/src/secondary_structure/pose.cpp
	/RNAMake/src/secondary_structure/residue.cc
	/RNAMake/src/secondary_structure/sequence_constraint.cc
	/RNAMake/src/secondary_structure/rna_structure.cpp
	/RNAMake/src/secondary_structure/secondary_structure_tree.cpp
	/RNAMake/src/secondary_structure/basepair.cc 
)
add_library(secondary_structure_lib STATIC ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib  -static)
###########################################################################
# eternabot
###########################################################################
set(eternabot_files
	/RNAMake/src/eternabot/feature_generator.cpp
	/RNAMake/src/eternabot/scorer.cpp
	/RNAMake/src/eternabot/sequence_designer.cpp
	/RNAMake/src/eternabot/strategy/berex_test.cpp
	/RNAMake/src/eternabot/strategy.cpp
	/RNAMake/src/eternabot/strategy/a_basic_test.cpp
	/RNAMake/src/eternabot/strategy/direction_of_gc.cpp
	/RNAMake/src/eternabot/strategy/clear_plot.cpp
	/RNAMake/src/eternabot/strategy/num_of_yellow.cpp 
)
add_library(eternabot_lib STATIC ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib secondary_structure_lib  -static)
###########################################################################
# structure
###########################################################################
set(structure_files
	/RNAMake/src/structure/structure.cc
	/RNAMake/src/structure/atom.cc
	/RNAMake/src/structure/chain.fwd.h
	/RNAMake/src/structure/residue_type_set.cc
	/RNAMake/src/structure/residue_type_set_manager.h
	/RNAMake/src/structure/beads.cc
	/RNAMake/src/structure/exceptions.cc
	/RNAMake/src/structure/chain.cc
	/RNAMake/src/structure/basepair.cc
	/RNAMake/src/structure/pdb_parser.cc
	/RNAMake/src/structure/close_chain.cc
	/RNAMake/src/structure/rna_structure.cpp
	/RNAMake/src/structure/basepair_state.cc
	/RNAMake/src/structure/residue_type.cc
	/RNAMake/src/structure/is_equal.cpp
	/RNAMake/src/structure/residue.cc 
)
add_library(structure_lib STATIC ${structure_files})
target_link_libraries(structure_lib util_lib  -static)
###########################################################################
# motif
###########################################################################
set(motif_files
	/RNAMake/src/motif/motif_state_aligner.cpp
	/RNAMake/src/motif/pose_factory.cc
	/RNAMake/src/motif/motif_to_secondary_structure.cc
	/RNAMake/src/motif/motif_ensemble.cpp
	/RNAMake/src/motif/motif_state_ensemble.cpp
	/RNAMake/src/motif/motif.cc
	/RNAMake/src/motif/motif_scorer.cc
	/RNAMake/src/motif/motif_state.cpp
	/RNAMake/src/motif/motif_factory.cc
	/RNAMake/src/motif/pose.cc 
)
add_library(motif_lib STATIC ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib  -static)
###########################################################################
# motif_tools
###########################################################################
set(motif_tools_files
	/RNAMake/src/motif_tools/segmenter.cpp 
)
add_library(motif_tools_lib STATIC ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib  -static)
###########################################################################
# resources
###########################################################################
set(resources_files
	/RNAMake/src/resources/added_motif_library.cc
	/RNAMake/src/resources/motif_ensemble_sqlite_connection.cpp
	/RNAMake/src/resources/motif_sqlite_connection.cc
	/RNAMake/src/resources/sqlite_library.cc
	/RNAMake/src/resources/motif_sqlite_library.cc
	/RNAMake/src/resources/resource_manager.cc
	/RNAMake/src/resources/motif_state_ensemble_sqlite_library.cpp
	/RNAMake/src/resources/motif_state_sqlite_library.cpp 
)
add_library(resources_lib STATIC ${resources_files})
target_link_libraries(resources_lib motif_lib  -static)
###########################################################################
# motif_data_structure
###########################################################################
set(motif_data_structure_files
	/RNAMake/src/motif_data_structure/motif_state_ensemble_tree.cpp
	/RNAMake/src/motif_data_structure/motif_state_node.cpp
	/RNAMake/src/motif_data_structure/motif_topology.cpp
	/RNAMake/src/motif_data_structure/motif_state_graph.cpp
	/RNAMake/src/motif_data_structure/motif_tree.cpp
	/RNAMake/src/motif_data_structure/motif_state_ensemble_graph.cc
	/RNAMake/src/motif_data_structure/motif_state_tree.fwd.h
	/RNAMake/src/motif_data_structure/motif_merger.cpp
	/RNAMake/src/motif_data_structure/motif_state_tree.cpp
	/RNAMake/src/motif_data_structure/motif_graph.cpp
	/RNAMake/src/motif_data_structure/motif_connection.cpp 
)
add_library(motif_data_structure_lib STATIC ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib  -static)
###########################################################################
# thermo_fluctuation
###########################################################################
set(thermo_fluctuation_files
	/RNAMake/src/thermo_fluctuation/graph/simulation.cc
	/RNAMake/src/thermo_fluctuation/graph/sterics.cc
	/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	/RNAMake/src/thermo_fluctuation/thermo_fluc_simulation.cpp
	/RNAMake/src/thermo_fluctuation/graph/scorer.cc
	/RNAMake/src/thermo_fluctuation/thermo_fluc_scorer.cpp
	/RNAMake/src/thermo_fluctuation/thermo_fluc_sampler.cpp
	/RNAMake/src/thermo_fluctuation/graph/logging.cc
	/RNAMake/src/thermo_fluctuation/thermo_fluc_relax.cc
	/RNAMake/src/thermo_fluctuation/graph/sampler.cc 
)
add_library(thermo_fluctuation_lib STATIC ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib  -static)
###########################################################################
# motif_search
###########################################################################
set(motif_search_files
	/RNAMake/src/motif_search/motif_state_monte_carlo.cc
	/RNAMake/src/motif_search/exhaustive/scorer.cc
	/RNAMake/src/motif_search/monte_carlo/scorer.cc
	/RNAMake/src/motif_search/path_finding/search.cc
	/RNAMake/src/motif_search/path_finding/scorer.cc
	/RNAMake/src/motif_search/path_finding/selector.cc
	/RNAMake/src/motif_search/solution_topology.cc
	/RNAMake/src/motif_search/exhaustive/motif_state_enumerator.cc
	/RNAMake/src/motif_search/exhaustive/search.cc
	/RNAMake/src/motif_search/solution_filter.cc
	/RNAMake/src/motif_search/monte_carlo/search.cc
	/RNAMake/src/motif_search/search.cc
	/RNAMake/src/motif_search/problem.cc
	/RNAMake/src/motif_search/path_finding/node.cc 
)
add_library(motif_search_lib STATIC ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib  -static)
###########################################################################
# sequence_optimization
###########################################################################
set(sequence_optimization_files
	/RNAMake/src/sequence_optimization/sequence_optimizer_3d.cpp
	/RNAMake/src/sequence_optimization/sequence_optimizer.cpp 
)
add_library(sequence_optimization_lib STATIC ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib eternabot_lib  -static)
add_library(all_lib STATIC /RNAMake/src/main.cpp)
target_link_libraries(all_lib motif_tools_lib thermo_fluctuation_lib motif_search_lib sequence_optimization_lib -static)
####################################################################################################
####################################################################################################
# Unittest Declarations 
####################################################################################################
###########################################################################
# base Tests 
###########################################################################
	add_executable(string_unittest /RNAMake/unittests/base/string_unittest.cpp)
	target_link_libraries(string_unittest base_lib -static )
	add_executable(cl_option_unittest /RNAMake/unittests/base/cl_option_unittest.cpp)
	target_link_libraries(cl_option_unittest base_lib -static )
	add_executable(application_unittest /RNAMake/unittests/base/application_unittest.cpp)
	target_link_libraries(application_unittest base_lib -static )
	add_executable(log_unittest /RNAMake/unittests/base/log_unittest.cc)
	target_link_libraries(log_unittest base_lib -static )
	add_executable(option_unittest /RNAMake/unittests/base/option_unittest.cpp)
	target_link_libraries(option_unittest base_lib -static )
###########################################################################
# math Tests 
###########################################################################
	add_executable(hashing_unittest /RNAMake/unittests/math/hashing_unittest.cc)
	target_link_libraries(hashing_unittest math_lib -static )
	add_executable(quaternion_unittest /RNAMake/unittests/math/quaternion_unittest.cc)
	target_link_libraries(quaternion_unittest math_lib -static )
	add_executable(rotation_unittest /RNAMake/unittests/math/rotation_unittest.cc)
	target_link_libraries(rotation_unittest math_lib -static )
	add_executable(xyz_matrix_unittest /RNAMake/unittests/math/xyz_matrix_unittest.cpp)
	target_link_libraries(xyz_matrix_unittest math_lib -static )
###########################################################################
# data_structure Tests 
###########################################################################
	add_executable(graph_unittest /RNAMake/unittests/data_structure/graph_unittest.cc)
	target_link_libraries(graph_unittest data_structure_lib -static )
	add_executable(tree_unittest /RNAMake/unittests/data_structure/tree_unittest.cpp)
	target_link_libraries(tree_unittest data_structure_lib -static )
	add_executable(graph_iter_list_unittest /RNAMake/unittests/data_structure/graph_iter_list_unittest.cc)
	target_link_libraries(graph_iter_list_unittest data_structure_lib -static )
	add_executable(graph_adjacency_list_unittest /RNAMake/unittests/data_structure/graph_adjacency_list_unittest.cc)
	target_link_libraries(graph_adjacency_list_unittest data_structure_lib -static )
###########################################################################
# util Tests 
###########################################################################
	add_executable(sqlite3_connection_unittest /RNAMake/unittests/util/sqlite3_connection_unittest.cpp)
	target_link_libraries(sqlite3_connection_unittest util_lib -static )
	add_executable(x3dna_unittest /RNAMake/unittests/util/x3dna_unittest.cpp)
	target_link_libraries(x3dna_unittest util_lib -static )
	add_executable(csv_unittest /RNAMake/unittests/util/csv_unittest.cpp)
	target_link_libraries(csv_unittest util_lib -static )
	add_executable(steric_lookup_unittest /RNAMake/unittests/util/steric_lookup_unittest.cpp)
	target_link_libraries(steric_lookup_unittest util_lib -static )
	add_executable(uuid_unittest /RNAMake/unittests/util/uuid_unittest.cpp)
	target_link_libraries(uuid_unittest util_lib -static )
###########################################################################
# vienna Tests 
###########################################################################
###########################################################################
# secondary_structure Tests 
###########################################################################
	add_executable(ss_structure_unittest /RNAMake/unittests/secondary_structure/ss_structure_unittest.cpp)
	target_link_libraries(ss_structure_unittest secondary_structure_lib -static )
	add_executable(secondary_structure_parser_unittest /RNAMake/unittests/secondary_structure/secondary_structure_parser_unittest.cpp)
	target_link_libraries(secondary_structure_parser_unittest secondary_structure_lib -static )
	add_executable(ss_motif_unittest /RNAMake/unittests/secondary_structure/ss_motif_unittest.cpp)
	target_link_libraries(ss_motif_unittest secondary_structure_lib -static )
	add_executable(ss_pose_unittest /RNAMake/unittests/secondary_structure/ss_pose_unittest.cpp)
	target_link_libraries(ss_pose_unittest secondary_structure_lib -static )
	add_executable(ss_chain_unittest /RNAMake/unittests/secondary_structure/ss_chain_unittest.cpp)
	target_link_libraries(ss_chain_unittest secondary_structure_lib -static )
	add_executable(ss_residue_unittest /RNAMake/unittests/secondary_structure/ss_residue_unittest.cpp)
	target_link_libraries(ss_residue_unittest secondary_structure_lib -static )
###########################################################################
# eternabot Tests 
###########################################################################
###########################################################################
# structure Tests 
###########################################################################
	add_executable(chain_closure_unittest /RNAMake/unittests/structure/chain_closure_unittest.cc)
	target_link_libraries(chain_closure_unittest structure_lib -static )
	add_executable(pdb_parser_unittest /RNAMake/unittests/structure/pdb_parser_unittest.cpp)
	target_link_libraries(pdb_parser_unittest structure_lib -static )
	add_executable(atom_unittest /RNAMake/unittests/structure/atom_unittest.cpp)
	target_link_libraries(atom_unittest structure_lib -static )
	add_executable(basepair_state_unittest /RNAMake/unittests/structure/basepair_state_unittest.cpp)
	target_link_libraries(basepair_state_unittest structure_lib -static )
	add_executable(residue_type_unittest /RNAMake/unittests/structure/residue_type_unittest.cpp)
	target_link_libraries(residue_type_unittest structure_lib -static )
	add_executable(basepair_unittest /RNAMake/unittests/structure/basepair_unittest.cpp)
	target_link_libraries(basepair_unittest structure_lib -static )
	add_executable(structure_unittest /RNAMake/unittests/structure/structure_unittest.cpp)
	target_link_libraries(structure_unittest structure_lib -static )
	add_executable(chain_unittest /RNAMake/unittests/structure/chain_unittest.cpp)
	target_link_libraries(chain_unittest structure_lib -static )
	add_executable(residue_unittest /RNAMake/unittests/structure/residue_unittest.cpp)
	target_link_libraries(residue_unittest structure_lib -static )
###########################################################################
# motif Tests 
###########################################################################
	add_executable(motif_state_unittest /RNAMake/unittests/motif/motif_state_unittest.cpp)
	target_link_libraries(motif_state_unittest motif_lib -static )
	add_executable(motif_to_secondary_structure_unittest /RNAMake/unittests/motif/motif_to_secondary_structure_unittest.cpp)
	target_link_libraries(motif_to_secondary_structure_unittest motif_lib -static )
	add_executable(motif_scorer_unittest /RNAMake/unittests/motif/motif_scorer_unittest.cpp)
	target_link_libraries(motif_scorer_unittest motif_lib -static )
	add_executable(motif_state_ensemble_unittest /RNAMake/unittests/motif/motif_state_ensemble_unittest.cpp)
	target_link_libraries(motif_state_ensemble_unittest motif_lib -static )
	add_executable(motif_unittest /RNAMake/unittests/motif/motif_unittest.cpp)
	target_link_libraries(motif_unittest motif_lib -static )
	add_executable(motif_factory_unittest /RNAMake/unittests/motif/motif_factory_unittest.cpp)
	target_link_libraries(motif_factory_unittest motif_lib -static )
	add_executable(motif_ensemble_unittest /RNAMake/unittests/motif/motif_ensemble_unittest.cpp)
	target_link_libraries(motif_ensemble_unittest motif_lib -static )
###########################################################################
# motif_tools Tests 
###########################################################################
	add_executable(segmenter_unittest /RNAMake/unittests/motif_tools/segmenter_unittest.cpp)
	target_link_libraries(segmenter_unittest motif_tools_lib -static )
###########################################################################
# resources Tests 
###########################################################################
	add_executable(motif_ensemble_sqlite_connection_unittest /RNAMake/unittests/resources/motif_ensemble_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_ensemble_sqlite_connection_unittest resources_lib -static )
	add_executable(resource_mananger_unittest /RNAMake/unittests/resources/resource_mananger_unittest.cpp)
	target_link_libraries(resource_mananger_unittest resources_lib -static )
	add_executable(motif_state_ensemble_sqlite_library_unittest /RNAMake/unittests/resources/motif_state_ensemble_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_ensemble_sqlite_library_unittest resources_lib -static )
	add_executable(added_motif_library_unittest /RNAMake/unittests/resources/added_motif_library_unittest.cpp)
	target_link_libraries(added_motif_library_unittest resources_lib -static )
	add_executable(motif_state_sqlite_library_unittest /RNAMake/unittests/resources/motif_state_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_sqlite_library_unittest resources_lib -static )
	add_executable(motif_sqlite_connection_unittest /RNAMake/unittests/resources/motif_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_sqlite_connection_unittest resources_lib -static )
	add_executable(motif_sqlite_library_unittest /RNAMake/unittests/resources/motif_sqlite_library_unittest.cpp)
	target_link_libraries(motif_sqlite_library_unittest resources_lib -static )
###########################################################################
# motif_data_structure Tests 
###########################################################################
	add_executable(motif_state_tree_unittest /RNAMake/unittests/motif_data_structure/motif_state_tree_unittest.cpp)
	target_link_libraries(motif_state_tree_unittest motif_data_structure_lib -static )
	add_executable(motif_graph_unittest /RNAMake/unittests/motif_data_structure/motif_graph_unittest.cpp)
	target_link_libraries(motif_graph_unittest motif_data_structure_lib -static )
	add_executable(motif_merger_unittest /RNAMake/unittests/motif_data_structure/motif_merger_unittest.cpp)
	target_link_libraries(motif_merger_unittest motif_data_structure_lib -static )
	add_executable(motif_state_ensemble_tree_unittest /RNAMake/unittests/motif_data_structure/motif_state_ensemble_tree_unittest.cpp)
	target_link_libraries(motif_state_ensemble_tree_unittest motif_data_structure_lib -static )
	add_executable(motif_topology_unittest /RNAMake/unittests/motif_data_structure/motif_topology_unittest.cpp)
	target_link_libraries(motif_topology_unittest motif_data_structure_lib -static )
	add_executable(motif_tree_unittest /RNAMake/unittests/motif_data_structure/motif_tree_unittest.cpp)
	target_link_libraries(motif_tree_unittest motif_data_structure_lib -static )
	add_executable(motif_state_graph_unittest /RNAMake/unittests/motif_data_structure/motif_state_graph_unittest.cpp)
	target_link_libraries(motif_state_graph_unittest motif_data_structure_lib -static )
###########################################################################
# thermo_fluctuation Tests 
###########################################################################
	add_executable(thermo_fluc_simulation_unittest /RNAMake/unittests/thermo_fluctuation/thermo_fluc_simulation_unittest.cpp)
	target_link_libraries(thermo_fluc_simulation_unittest thermo_fluctuation_lib -static )
	add_executable(thermo_fluc_sampler_unittest /RNAMake/unittests/thermo_fluctuation/thermo_fluc_sampler_unittest.cpp)
	target_link_libraries(thermo_fluc_sampler_unittest thermo_fluctuation_lib -static )
###########################################################################
# motif_search Tests 
###########################################################################
	add_executable(monte_carlo_search_unittest /RNAMake/unittests/motif_search/monte_carlo_search_unittest.cc)
	target_link_libraries(monte_carlo_search_unittest motif_search_lib -static )
	add_executable(exhaustive_search_unittest /RNAMake/unittests/motif_search/exhaustive_search_unittest.cc)
	target_link_libraries(exhaustive_search_unittest motif_search_lib -static )
	add_executable(motif_state_search_unittest /RNAMake/unittests/motif_search/motif_state_search_unittest.cpp)
	target_link_libraries(motif_state_search_unittest motif_search_lib -static )
###########################################################################
# sequence_optimization Tests 
###########################################################################
	add_executable(sequence_optimizer_unittest /RNAMake/unittests/sequence_optimization/sequence_optimizer_unittest.cpp)
	target_link_libraries(sequence_optimizer_unittest sequence_optimization_lib -static )
####################################################################################################
####################################################################################################
# App declarations
###############################################################################################################################################################################
# ddg_tecto
###########################################################################
	add_executable(ddg_tecto /RNAMake/apps/ddg_tecto/ddg_tecto.cc)
	target_link_libraries(ddg_tecto all_lib -static )
###########################################################################
# simulate_tectos_devel
###########################################################################
	add_executable(simulate_tectos_devel /RNAMake/apps/simulate_tectos_devel/simulate_tectos_devel.cc)
	target_link_libraries(simulate_tectos_devel all_lib -static )
###########################################################################
# sequence_optimizer_app
###########################################################################
	add_executable(sequence_optimizer_app /RNAMake/apps/sequence_optimizer_app/sequence_optimizer_app.cpp)
	target_link_libraries(sequence_optimizer_app all_lib -static )
###########################################################################
# thermo_simulation
###########################################################################
	add_executable(thermo_simulation /RNAMake/apps/thermo_simulation/thermo_simulation.cpp)
	target_link_libraries(thermo_simulation all_lib -static )
###########################################################################
# sequence_optimization_benchmarks
###########################################################################
	add_executable(sequence_optimization_benchmarks /RNAMake/apps/sequence_optimization_benchmarks/sequence_optimization_benchmarks.cc)
	target_link_libraries(sequence_optimization_benchmarks all_lib -static )
###########################################################################
# design_rna_scaffold
###########################################################################
	add_executable(design_rna_scaffold /RNAMake/apps/design_rna_scaffold/design_rna_scaffold.cc)
	target_link_libraries(design_rna_scaffold all_lib -static )
###########################################################################
# get_best_solutions
###########################################################################
	add_executable(get_best_solutions /RNAMake/apps/get_best_solutions/get_best_solutions.cc)
	target_link_libraries(get_best_solutions all_lib -static )
###########################################################################
# eternabot
###########################################################################
	add_executable(eternabot /RNAMake/apps/eternabot/eternabot.cpp)
	target_link_libraries(eternabot all_lib -static )
####################################################################################################
