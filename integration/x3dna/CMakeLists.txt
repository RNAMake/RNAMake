####################################################################################################
# Project Level Info
####################################################################################################
cmake_minimum_required(VERSION 3.0)
set(CMAKE_BUILD_TYPE Release)
project(RNAMake)
set(CMAKE_CXX_STANDARD 14)
set(BUILD DEBUG)

set(RNAMAKE /Users/cjurich/projects/RNAMake/)

include(${RNAMAKE}/cmake/build/compiler.cmake)
include(${RNAMAKE}/cmake/build/sqlite.cmake)
include_directories(${RNAMAKE}/src/external/)
include_directories(${RNAMAKE}/unittests/)
include_directories(${RNAMAKE}/apps/)
include_directories(${RNAMAKE}/src/)
find_library(base_lib . )
#add_subdirectory(${RNAMAKE}/src/)
#add_subdirectory(${RNAMAKE}/src/)
#include(${RNAMAKE}/src/base/CMakeLists.txt)
#find_package(base_lib
#        PATHS
#        ../../src/base/
#    )
####################################################################################################
####################################################################################################
# Library declarations
####################################################################################################
###########################################################################
# base
###########################################################################
set(base_files
	${RNAMAKE}/src/base/string.cc
	${RNAMAKE}/src/base/application.cpp
	${RNAMAKE}/src/base/file_io.cc
	${RNAMAKE}/src/base/cl_option.cc
	${RNAMAKE}/src/base/assertions.h
	${RNAMAKE}/src/base/settings.cc
	${RNAMAKE}/src/base/vector_container.h
	${RNAMAKE}/src/base/option.cc
	${RNAMAKE}/src/base/util.hpp
	${RNAMAKE}/src/base/sys_interface.cpp
	${RNAMAKE}/src/base/backtrace.cpp
	${RNAMAKE}/src/base/command_line_parser.cpp
	${RNAMAKE}/src/base/global_constants.cc
	${RNAMAKE}/src/base/types.h
	${RNAMAKE}/src/base/log.cc
)
add_library(base_lib  ${base_files})
target_link_libraries(base_lib  )
###########################################################################
# math
###########################################################################
set(math_files
	${RNAMAKE}/src/math/xyz_matrix.h
	${RNAMAKE}/src/math/euler.cc
	${RNAMAKE}/src/math/xyz_vector.h
	${RNAMAKE}/src/math/quaternion.cc
	${RNAMAKE}/src/math/stats.cc
	${RNAMAKE}/src/math/numerical.cc
	${RNAMAKE}/src/math/transform.h
	${RNAMAKE}/src/math/hashing.cc
)
add_library(math_lib  ${math_files})
target_link_libraries(math_lib base_lib  )
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	${RNAMAKE}/src/data_structure/tree/tree.cc
	${RNAMAKE}/src/data_structure/graph/graph_node.fwd.h
	${RNAMAKE}/src/data_structure/graph_base.cc
	${RNAMAKE}/src/data_structure/graph/graph.cc
	${RNAMAKE}/src/data_structure/tree/tree_node.cc
	${RNAMAKE}/src/data_structure/graph_adjacency_list.cc
	${RNAMAKE}/src/data_structure/graph_iter_list.cc
	${RNAMAKE}/src/data_structure/graph/graph_node.cc
	${RNAMAKE}/src/data_structure/graph.cc
)
add_library(data_structure_lib  ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  )
###########################################################################
# util
###########################################################################
set(util_files
	${RNAMAKE}/src/util/dssr.cpp
	${RNAMAKE}/src/util/monte_carlo.cpp
	${RNAMAKE}/src/util/basic_io.cpp
	${RNAMAKE}/src/util/random_number_generator.h
	${RNAMAKE}/src/util/csv.cpp
	${RNAMAKE}/src/util/steric_lookup.cpp
	${RNAMAKE}/src/util/x3dna.cc
	${RNAMAKE}/src/util/sqlite3_connection.cc
	${RNAMAKE}/src/util/motif_type.cc
	${RNAMAKE}/src/util/uuid.cc
	${RNAMAKE}/src/util/cartesian_product.cc
)
add_library(util_lib  ${util_files})
target_link_libraries(util_lib math_lib sqlite3 )
###########################################################################
# vienna
###########################################################################
set(vienna_files
	${RNAMAKE}/src/vienna/pair_mat.h
	${RNAMAKE}/src/vienna/energy_par.cc
	${RNAMAKE}/src/vienna/vienna.cc
)
add_library(vienna_lib  ${vienna_files})
target_link_libraries(vienna_lib base_lib  )
###########################################################################
# secondary_structure
###########################################################################
set(secondary_structure_files
	${RNAMAKE}/src/secondary_structure/residue.cc
	${RNAMAKE}/src/secondary_structure/motif.cc
	${RNAMAKE}/src/secondary_structure/util.cpp
	${RNAMAKE}/src/secondary_structure/sequence_constraint.cc
	${RNAMAKE}/src/secondary_structure/secondary_structure_tree.cpp
	${RNAMAKE}/src/secondary_structure/sequence_tools.cc
	${RNAMAKE}/src/secondary_structure/basepair.cc
	${RNAMAKE}/src/secondary_structure/rna_structure.cpp
	${RNAMAKE}/src/secondary_structure/structure.cpp
	${RNAMAKE}/src/secondary_structure/chain.cc
	${RNAMAKE}/src/secondary_structure/secondary_structure_parser.cpp
	${RNAMAKE}/src/secondary_structure/pose.cpp
)
add_library(secondary_structure_lib  ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib  )
###########################################################################
# eternabot
###########################################################################
set(eternabot_files
	${RNAMAKE}/src/eternabot/strategy/a_basic_test.cpp
	${RNAMAKE}/src/eternabot/strategy.cpp
	${RNAMAKE}/src/eternabot/scorer.cpp
	${RNAMAKE}/src/eternabot/strategy/clear_plot.cpp
	${RNAMAKE}/src/eternabot/strategy/direction_of_gc.cpp
	${RNAMAKE}/src/eternabot/strategy/berex_test.cpp
	${RNAMAKE}/src/eternabot/feature_generator.cpp
	${RNAMAKE}/src/eternabot/sequence_designer.cpp
	${RNAMAKE}/src/eternabot/strategy/num_of_yellow.cpp
)
add_library(eternabot_lib  ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib secondary_structure_lib  )
###########################################################################
# structure
###########################################################################
set(structure_files
	${RNAMAKE}/src/structure/residue_type_set.cc
	${RNAMAKE}/src/structure/close_chain.cc
	${RNAMAKE}/src/structure/exceptions.cc
	${RNAMAKE}/src/structure/residue_type.cc
	${RNAMAKE}/src/structure/atom.cc
	${RNAMAKE}/src/structure/structure.cc
	${RNAMAKE}/src/structure/cif_parser.cc
	${RNAMAKE}/src/structure/basepair_state.cc
	${RNAMAKE}/src/structure/chain.fwd.h
	${RNAMAKE}/src/structure/beads.cc
	${RNAMAKE}/src/structure/basepair.cc
	${RNAMAKE}/src/structure/chain.cc
	${RNAMAKE}/src/structure/residue.cc
	${RNAMAKE}/src/structure/is_equal.cpp
	${RNAMAKE}/src/structure/rna_structure.cpp
	${RNAMAKE}/src/structure/pdb_parser.cc
	${RNAMAKE}/src/structure/residue_type_set_manager.h
)
add_library(structure_lib  ${structure_files})
target_link_libraries(structure_lib util_lib  )
###########################################################################
# motif
###########################################################################
set(motif_files
	${RNAMAKE}/src/motif/motif_state_aligner.cpp
	${RNAMAKE}/src/motif/motif.cc
	${RNAMAKE}/src/motif/motif_scorer.cc
	${RNAMAKE}/src/motif/pose.cc
	${RNAMAKE}/src/motif/motif_state.cpp
	${RNAMAKE}/src/motif/motif_state_ensemble.cpp
	${RNAMAKE}/src/motif/pose_factory.cc
	${RNAMAKE}/src/motif/motif_ensemble.cpp
	${RNAMAKE}/src/motif/motif_factory.cc
	${RNAMAKE}/src/motif/motif_to_secondary_structure.cc
)
add_library(motif_lib  ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib  )
###########################################################################
# motif_tools
###########################################################################
set(motif_tools_files
	${RNAMAKE}/src/motif_tools/segmenter.cpp
)
add_library(motif_tools_lib  ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib  )
###########################################################################
# resources
###########################################################################
set(resources_files
	${RNAMAKE}/src/resources/motif_ensemble_sqlite_connection.cpp
	${RNAMAKE}/src/resources/motif_state_ensemble_sqlite_library.cpp
	${RNAMAKE}/src/resources/added_motif_library.cc
	${RNAMAKE}/src/resources/motif_sqlite_connection.cc
	${RNAMAKE}/src/resources/motif_sqlite_library.cc
	${RNAMAKE}/src/resources/motif_state_sqlite_library.cpp
	${RNAMAKE}/src/resources/resource_manager.cc
	${RNAMAKE}/src/resources/sqlite_library.cc
)
add_library(resources_lib  ${resources_files})
target_link_libraries(resources_lib motif_lib  )
###########################################################################
# motif_data_structure
###########################################################################
set(motif_data_structure_files
	${RNAMAKE}/src/motif_data_structure/motif_state_ensemble_graph.cc
	${RNAMAKE}/src/motif_data_structure/motif_state_tree.fwd.h
	${RNAMAKE}/src/motif_data_structure/motif_connection.cpp
	${RNAMAKE}/src/motif_data_structure/motif_state_tree.cpp
	${RNAMAKE}/src/motif_data_structure/motif_state_node.cpp
	${RNAMAKE}/src/motif_data_structure/motif_state_graph.cpp
	${RNAMAKE}/src/motif_data_structure/motif_graph.cpp
	${RNAMAKE}/src/motif_data_structure/motif_state_ensemble_tree.cpp
	${RNAMAKE}/src/motif_data_structure/motif_topology.cpp
	${RNAMAKE}/src/motif_data_structure/motif_tree.cpp
	${RNAMAKE}/src/motif_data_structure/motif_merger.cpp
)
add_library(motif_data_structure_lib  ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib  )
###########################################################################
# thermo_fluctuation
###########################################################################
set(thermo_fluctuation_files
	${RNAMAKE}/src/thermo_fluctuation/graph/sterics.cc
	${RNAMAKE}/src/thermo_fluctuation/thermo_fluc_sampler.cpp
	${RNAMAKE}/src/thermo_fluctuation/thermo_fluc_simulation.cpp
	${RNAMAKE}/src/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	${RNAMAKE}/src/thermo_fluctuation/graph/logging.cc
	${RNAMAKE}/src/thermo_fluctuation/thermo_fluc_relax.cc
	${RNAMAKE}/src/thermo_fluctuation/thermo_fluc_scorer.cpp
	${RNAMAKE}/src/thermo_fluctuation/graph/simulation.cc
	${RNAMAKE}/src/thermo_fluctuation/graph/scorer.cc
	${RNAMAKE}/src/thermo_fluctuation/graph/sampler.cc
)
add_library(thermo_fluctuation_lib  ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib  )
###########################################################################
# motif_search
###########################################################################
set(motif_search_files
	${RNAMAKE}/src/motif_search/exhaustive/motif_state_enumerator.cc
	${RNAMAKE}/src/motif_search/motif_state_monte_carlo.cc
	${RNAMAKE}/src/motif_search/path_finding/node.cc
	${RNAMAKE}/src/motif_search/path_finding/selector.cc
	${RNAMAKE}/src/motif_search/path_finding/scorer.cc
	${RNAMAKE}/src/motif_search/path_finding/search.cc
	${RNAMAKE}/src/motif_search/exhaustive/search.cc
	${RNAMAKE}/src/motif_search/monte_carlo/scorer.cc
	${RNAMAKE}/src/motif_search/problem.cc
	${RNAMAKE}/src/motif_search/monte_carlo/search.cc
	${RNAMAKE}/src/motif_search/search.cc
	${RNAMAKE}/src/motif_search/exhaustive/scorer.cc
	${RNAMAKE}/src/motif_search/solution_topology.cc
	${RNAMAKE}/src/motif_search/solution_filter.cc
)
add_library(motif_search_lib  ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib  )
###################################################################################################
# x3dna integration test
###################################################################################################

add_executable(x3dna_integration x3dna_integration.cpp)
target_link_libraries(x3dna_integration util_lib)
