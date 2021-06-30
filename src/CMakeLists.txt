###########################################################################
# base
###########################################################################

set(RNAMAKE_SRC_DIR ${RNAMAKE}/src/)
set(RNAMAKE_EXT_SRC_DIR ${RNAMAKE_SRC_DIR}/external/)

include_directories( ${RNAMAKE_SRC_DIR} )
include_directories( ${RNAMAKE_EXT_SRC_DIR} )

if ( NOT EXISTS ${BUILD_DIR}/lib )
    make_directory( ${BUILD_DIR}/lib )
endif()

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BUILD_DIR}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}/lib )

# building sqlite
include(${RNAMAKE_EXT_SRC_DIR}/sqlite/CMakeLists.txt)

set(base_files
	${RNAMAKE_SRC_DIR}/base/exception.cpp
	${RNAMAKE_SRC_DIR}/base/assertions.h
	${RNAMAKE_SRC_DIR}/base/file_io.cc
	${RNAMAKE_SRC_DIR}/base/sys_interface.cpp
	${RNAMAKE_SRC_DIR}/base/cl_option.cc
	${RNAMAKE_SRC_DIR}/base/option.cc
	${RNAMAKE_SRC_DIR}/base/global_constants.cc
	${RNAMAKE_SRC_DIR}/base/types.h
	${RNAMAKE_SRC_DIR}/base/backtrace.cpp
	${RNAMAKE_SRC_DIR}/base/command_line_parser.cpp
	${RNAMAKE_SRC_DIR}/base/string.cc
	${RNAMAKE_SRC_DIR}/base/vector_container.h
	${RNAMAKE_SRC_DIR}/base/log.cc
	${RNAMAKE_SRC_DIR}/base/util.hpp
	${RNAMAKE_SRC_DIR}/base/settings.cc
	${RNAMAKE_SRC_DIR}/base/env_manager.cpp
	${RNAMAKE_SRC_DIR}/base/application.cpp 
)
add_library(base_lib  ${base_files})
target_link_libraries(base_lib  )
###########################################################################
# math
###########################################################################
set(math_files
	${RNAMAKE_SRC_DIR}/math/hashing.cc
	${RNAMAKE_SRC_DIR}/math/xyz_matrix.h
	${RNAMAKE_SRC_DIR}/math/quaternion.cc
	${RNAMAKE_SRC_DIR}/math/stats.cc
	${RNAMAKE_SRC_DIR}/math/numerical.cc
	${RNAMAKE_SRC_DIR}/math/xyz_vector.h
	${RNAMAKE_SRC_DIR}/math/transform.h
	${RNAMAKE_SRC_DIR}/math/euler.cc 
)
add_library(math_lib  ${math_files})
target_link_libraries(math_lib base_lib  )
###########################################################################
# data_structure
###########################################################################
set(data_structure_files
	${RNAMAKE_SRC_DIR}/data_structure/tree/tree.cc
	${RNAMAKE_SRC_DIR}/data_structure/graph_adjacency_list.cc
	${RNAMAKE_SRC_DIR}/data_structure/graph_iter_list.cc
	${RNAMAKE_SRC_DIR}/data_structure/graph_base.cc
	${RNAMAKE_SRC_DIR}/data_structure/graph.cc
	${RNAMAKE_SRC_DIR}/data_structure/tree/tree_node.fwd.hh
	${RNAMAKE_SRC_DIR}/data_structure/graph/graph.cc
	${RNAMAKE_SRC_DIR}/data_structure/graph/graph_node.cc
	${RNAMAKE_SRC_DIR}/data_structure/tree/tree_node.cc 
)
add_library(data_structure_lib  ${data_structure_files})
target_link_libraries(data_structure_lib base_lib  )
###########################################################################
# util
###########################################################################
set(util_files
	${RNAMAKE_SRC_DIR}/util/nrutil.cpp
	${RNAMAKE_SRC_DIR}/util/string_util.h
	${RNAMAKE_SRC_DIR}/util/app_fncs.cpp
	${RNAMAKE_SRC_DIR}/util/cartesian_product.cc
	${RNAMAKE_SRC_DIR}/util/sqlite3_connection.cc
	${RNAMAKE_SRC_DIR}/util/monte_carlo.cpp
	${RNAMAKE_SRC_DIR}/util/steric_lookup.cpp
	${RNAMAKE_SRC_DIR}/util/visit_struct.h
	${RNAMAKE_SRC_DIR}/util/motif_type.cc
	${RNAMAKE_SRC_DIR}/util/uuid.cc
	${RNAMAKE_SRC_DIR}/util/dssr.cpp
	${RNAMAKE_SRC_DIR}/util/basic_io.cpp
	${RNAMAKE_SRC_DIR}/util/x3dna_fncs.h
	${RNAMAKE_SRC_DIR}/util/find_pair.cpp
	${RNAMAKE_SRC_DIR}/util/ana_fncs.cpp
	${RNAMAKE_SRC_DIR}/util/fncs_slre.cpp
	${RNAMAKE_SRC_DIR}/util/testing_tools.cpp
	${RNAMAKE_SRC_DIR}/util/visit_struct_intrusive.h
	${RNAMAKE_SRC_DIR}/util/x3dna.cc
	${RNAMAKE_SRC_DIR}/util/csv.h
	${RNAMAKE_SRC_DIR}/util/cmn_fncs.cpp
	${RNAMAKE_SRC_DIR}/util/analyze.cpp
	${RNAMAKE_SRC_DIR}/util/random_number_generator.h
	${RNAMAKE_SRC_DIR}/util/x3dna_src.h 
)
add_library(util_lib  ${util_files})
target_link_libraries(util_lib math_lib sqlite3 )
###########################################################################
# vienna
###########################################################################
set(vienna_files
	${RNAMAKE_SRC_DIR}/vienna/pair_mat.h
	${RNAMAKE_SRC_DIR}/vienna/vienna.cc
	${RNAMAKE_SRC_DIR}/vienna/energy_par.cc 
)
add_library(vienna_lib  ${vienna_files})
target_link_libraries(vienna_lib base_lib  )
###########################################################################
# secondary_structure
###########################################################################
set(secondary_structure_files
	${RNAMAKE_SRC_DIR}/secondary_structure/util.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/secondary_structure_tree.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/sequence_constraint.cc
	${RNAMAKE_SRC_DIR}/secondary_structure/motif.cc
	${RNAMAKE_SRC_DIR}/secondary_structure/residue.cc
	${RNAMAKE_SRC_DIR}/secondary_structure/chain.cc
	${RNAMAKE_SRC_DIR}/secondary_structure/basepair.cc
	${RNAMAKE_SRC_DIR}/secondary_structure/rna_structure.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/pose.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/secondary_structure_parser.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/structure.cpp
	${RNAMAKE_SRC_DIR}/secondary_structure/sequence_tools.cc 
)
add_library(secondary_structure_lib  ${secondary_structure_files})
target_link_libraries(secondary_structure_lib util_lib  )
###########################################################################
# eternabot
###########################################################################
set(eternabot_files
	${RNAMAKE_SRC_DIR}/eternabot/sequence_designer.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy/clear_plot.cpp
	${RNAMAKE_SRC_DIR}/eternabot/feature_generator.cpp
	${RNAMAKE_SRC_DIR}/eternabot/scorer.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy/direction_of_gc.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy/berex_test.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy/num_of_yellow.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy.cpp
	${RNAMAKE_SRC_DIR}/eternabot/strategy/a_basic_test.cpp 
)
add_library(eternabot_lib  ${eternabot_files})
target_link_libraries(eternabot_lib vienna_lib secondary_structure_lib  )
###########################################################################
# structure
###########################################################################
set(structure_files
	${RNAMAKE_SRC_DIR}/structure/residue_type_set_manager.h
	${RNAMAKE_SRC_DIR}/structure/basepair.cc
	${RNAMAKE_SRC_DIR}/structure/cif_parser.cc
	${RNAMAKE_SRC_DIR}/structure/is_equal.cpp
	${RNAMAKE_SRC_DIR}/structure/basepair_state.cc
	${RNAMAKE_SRC_DIR}/structure/close_chain.cc
	${RNAMAKE_SRC_DIR}/structure/basepair_state.fwd.h
	${RNAMAKE_SRC_DIR}/structure/beads.cc
	${RNAMAKE_SRC_DIR}/structure/residue.cc
	${RNAMAKE_SRC_DIR}/structure/pdb_parser.cc
	${RNAMAKE_SRC_DIR}/structure/exceptions.cc
	${RNAMAKE_SRC_DIR}/structure/residue_type.cc
	${RNAMAKE_SRC_DIR}/structure/structure.cc
	${RNAMAKE_SRC_DIR}/structure/residue_type_set.cc
	${RNAMAKE_SRC_DIR}/structure/rna_structure.cpp
	${RNAMAKE_SRC_DIR}/structure/atom.cc
	${RNAMAKE_SRC_DIR}/structure/chain.cc 
)
add_library(structure_lib  ${structure_files})
target_link_libraries(structure_lib util_lib  )
###########################################################################
# motif
###########################################################################
set(motif_files
	${RNAMAKE_SRC_DIR}/motif/motif_state.cpp
	${RNAMAKE_SRC_DIR}/motif/motif_factory.cc
	${RNAMAKE_SRC_DIR}/motif/motif_to_secondary_structure.cc
	${RNAMAKE_SRC_DIR}/motif/motif_scorer.cc
	${RNAMAKE_SRC_DIR}/motif/motif_state_aligner.cpp
	${RNAMAKE_SRC_DIR}/motif/motif_ensemble.cpp
	${RNAMAKE_SRC_DIR}/motif/pose_factory.cc
	${RNAMAKE_SRC_DIR}/motif/motif_state_ensemble.cpp
	${RNAMAKE_SRC_DIR}/motif/motif.cc
	${RNAMAKE_SRC_DIR}/motif/pose.cc 
)
add_library(motif_lib  ${motif_files})
target_link_libraries(motif_lib structure_lib secondary_structure_lib  )
###########################################################################
# motif_tools
###########################################################################
set(motif_tools_files
	${RNAMAKE_SRC_DIR}/motif_tools/segmenter.cpp 
)
add_library(motif_tools_lib  ${motif_tools_files})
target_link_libraries(motif_tools_lib motif_lib  )
###########################################################################
# resources
###########################################################################
set(resources_files
	${RNAMAKE_SRC_DIR}/resources/resource_manager.cc
	${RNAMAKE_SRC_DIR}/resources/motif_state_sqlite_library.cpp
	${RNAMAKE_SRC_DIR}/resources/motif_ensemble_sqlite_connection.cpp
	${RNAMAKE_SRC_DIR}/resources/added_motif_library.cc
	${RNAMAKE_SRC_DIR}/resources/sqlite_library.cc
	${RNAMAKE_SRC_DIR}/resources/motif_state_ensemble_sqlite_library.cpp
	${RNAMAKE_SRC_DIR}/resources/motif_sqlite_library.cc
	${RNAMAKE_SRC_DIR}/resources/motif_sqlite_connection.cc 
)
add_library(resources_lib  ${resources_files})
target_link_libraries(resources_lib motif_lib  )
###########################################################################
# motif_data_structure
###########################################################################
set(motif_data_structure_files
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_ensemble_tree.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_topology.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_tree.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_tree.fwd.h
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_node.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_tree.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_graph.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_graph.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_state_ensemble_graph.cc
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_merger.cpp
	${RNAMAKE_SRC_DIR}/motif_data_structure/motif_connection.cpp 
)
add_library(motif_data_structure_lib  ${motif_data_structure_files})
target_link_libraries(motif_data_structure_lib resources_lib data_structure_lib  )
###########################################################################
# thermo_fluctuation
###########################################################################
set(thermo_fluctuation_files
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/graph/sterics.cc
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/thermo_fluc_scorer.cpp
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/graph/scorer.cc
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/graph/sampler.cc
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/graph/logging.cc
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/thermo_fluc_simulation.cpp
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/thermo_fluc_relax.cc
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/thermo_fluc_simulation_devel.cpp
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/thermo_fluc_sampler.cpp
	${RNAMAKE_SRC_DIR}/thermo_fluctuation/graph/simulation.cc 
)
add_library(thermo_fluctuation_lib  ${thermo_fluctuation_files})
target_link_libraries(thermo_fluctuation_lib motif_data_structure_lib  )
###########################################################################
# motif_search
###########################################################################
set(motif_search_files
	${RNAMAKE_SRC_DIR}/motif_search/path_finding/search.cc
	${RNAMAKE_SRC_DIR}/motif_search/path_finding/node.cc
	${RNAMAKE_SRC_DIR}/motif_search/path_finding/selector.cc
	${RNAMAKE_SRC_DIR}/motif_search/problem.cc
	${RNAMAKE_SRC_DIR}/motif_search/motif_state_monte_carlo.cc
	${RNAMAKE_SRC_DIR}/motif_search/monte_carlo/search.cc
	${RNAMAKE_SRC_DIR}/motif_search/search.cc
	${RNAMAKE_SRC_DIR}/motif_search/solution_topology.cc
	${RNAMAKE_SRC_DIR}/motif_search/exhaustive/search.cc
	${RNAMAKE_SRC_DIR}/motif_search/solution_filter.cc
	${RNAMAKE_SRC_DIR}/motif_search/monte_carlo/scorer.cc
	${RNAMAKE_SRC_DIR}/motif_search/exhaustive/scorer.cc
	${RNAMAKE_SRC_DIR}/motif_search/exhaustive/motif_state_enumerator.cc
	${RNAMAKE_SRC_DIR}/motif_search/path_finding/scorer.cc 
)
add_library(motif_search_lib  ${motif_search_files})
target_link_libraries(motif_search_lib motif_data_structure_lib  )
###########################################################################
# sequence_optimization
###########################################################################
set(sequence_optimization_files
	${RNAMAKE_SRC_DIR}/sequence_optimization/sequence_optimizer.cpp
	${RNAMAKE_SRC_DIR}/sequence_optimization/sequence_optimizer_3d.cpp 
)
add_library(sequence_optimization_lib  ${sequence_optimization_files})
target_link_libraries(sequence_optimization_lib motif_data_structure_lib eternabot_lib  )

add_library(RNAMake  ${RNAMAKE_SRC_DIR}/main.cpp)
target_link_libraries( RNAMake motif_tools_lib thermo_fluctuation_lib motif_search_lib sequence_optimization_lib )

