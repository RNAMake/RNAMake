set(RNAMAKE_UNITTEST_DIR ${RNAMAKE}/unittests/ )

include_directories( ${RNAMAKE_UNITTEST_DIR} )

if ( NOT EXISTS ${BUILD_DIR}/unittests )
    make_directory( ${BUILD_DIR}/unittests )
endif()

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BUILD_DIR}/unittests )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD_DIR}/unittests )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}/unittests )


###########################################################################
# base Tests 
###########################################################################
	add_executable(env_manager_unittest ${RNAMAKE_UNITTEST_DIR}/base/env_manager_unittest.cpp)
	target_link_libraries(env_manager_unittest base_lib  )
    add_test(env_manager_unittest env_manager_unittest)
	
    add_executable(application_unittest ${RNAMAKE_UNITTEST_DIR}/base/application_unittest.cpp)
    target_link_libraries(application_unittest base_lib  )
	add_test(application_unittest application_unittest)
	
    add_executable(cl_option_unittest ${RNAMAKE_UNITTEST_DIR}/base/cl_option_unittest.cpp)
	target_link_libraries(cl_option_unittest base_lib  )
	add_test(cl_option_unittest cl_option_unittest)

    add_executable(option_unittest ${RNAMAKE_UNITTEST_DIR}/base/option_unittest.cpp)
	target_link_libraries(option_unittest base_lib  )
	add_test(option_unittest option_unittest)
	
    add_executable(string_unittest ${RNAMAKE_UNITTEST_DIR}/base/string_unittest.cpp)
	target_link_libraries(string_unittest base_lib  )
	add_test(string_unittest string_unittest)
	
    add_executable(log_unittest ${RNAMAKE_UNITTEST_DIR}/base/log_unittest.cpp)
	target_link_libraries(log_unittest base_lib  )
	add_test(log_unittest log_unittest)
###########################################################################
# math Tests 
###########################################################################
	add_executable(hashing_unittest ${RNAMAKE_UNITTEST_DIR}/math/hashing_unittest.cpp)
	target_link_libraries(hashing_unittest math_lib  )
	add_test(hashing_unittest hashing_unittest)
	
    add_executable(quaternion_unittest ${RNAMAKE_UNITTEST_DIR}/math/quaternion_unittest.cpp)
	target_link_libraries(quaternion_unittest math_lib  )
	add_test(quaternion_unittest quaternion_unittest)
	
    add_executable(xyz_matrix_unittest ${RNAMAKE_UNITTEST_DIR}/math/xyz_matrix_unittest.cpp)
	target_link_libraries(xyz_matrix_unittest math_lib  )
	add_test(xyz_matrix_unittest xyz_matrix_unittest)
	
    add_executable(rotation_unittest ${RNAMAKE_UNITTEST_DIR}/math/rotation_unittest.cpp)
	target_link_libraries(rotation_unittest math_lib  )
	add_test(rotation_unittest rotation_unittest)
###########################################################################
# data_structure Tests 
###########################################################################
	add_executable(tree_unittest ${RNAMAKE_UNITTEST_DIR}/data_structure/tree_unittest.cpp)
	target_link_libraries(tree_unittest data_structure_lib  )
	add_test(tree_unittest tree_unittest)
	
    add_executable(graph_iter_list_unittest ${RNAMAKE_UNITTEST_DIR}/data_structure/graph_iter_list_unittest.cpp)
	target_link_libraries(graph_iter_list_unittest data_structure_lib  )
	add_test(graph_iter_list_unittest graph_iter_list_unittest)
	
    add_executable(graph_adjacency_list_unittest ${RNAMAKE_UNITTEST_DIR}/data_structure/graph_adjacency_list_unittest.cpp)
	target_link_libraries(graph_adjacency_list_unittest data_structure_lib  )
	add_test(graph_adjacency_list_unittest graph_adjacency_list_unittest)
	
    add_executable(graph_unittest ${RNAMAKE_UNITTEST_DIR}/data_structure/graph_unittest.cpp)
	target_link_libraries(graph_unittest data_structure_lib  ) 
    add_test(graph_unittest graph_unittest) 
    ###########################################################################
# util Tests 
###########################################################################
	add_executable(steric_lookup_unittest ${RNAMAKE_UNITTEST_DIR}/util/steric_lookup_unittest.cpp)
	target_link_libraries(steric_lookup_unittest util_lib  )
	add_test(steric_lookup_unittest steric_lookup_unittest)
	
    add_executable(csv_unittest ${RNAMAKE_UNITTEST_DIR}/util/csv_unittest.cpp)
	target_link_libraries(csv_unittest util_lib  )
	add_test(csv_unittest csv_unittest)
	
    add_executable(sqlite3_connection_unittest ${RNAMAKE_UNITTEST_DIR}/util/sqlite3_connection_unittest.cpp)
	target_link_libraries(sqlite3_connection_unittest util_lib  )
	add_test(sqlite3_connection_unittest sqlite3_connection_unittest)
	
    add_executable(dssr_unittest ${RNAMAKE_UNITTEST_DIR}/util/dssr_unittest.cpp)
	target_link_libraries(dssr_unittest util_lib  )
	add_test(dssr_unittest dssr_unittest)
	
    add_executable(uuid_unittest ${RNAMAKE_UNITTEST_DIR}/util/uuid_unittest.cpp)
	target_link_libraries(uuid_unittest util_lib  )
	add_test(uuid_unittest uuid_unittest)
	
    add_executable(motif_type_unittest ${RNAMAKE_UNITTEST_DIR}/util/motif_type_unittest.cpp)
	target_link_libraries(motif_type_unittest util_lib  )
	add_test(motif_type_unittest motif_type_unittest)
###########################################################################
# vienna Tests 
###########################################################################
###########################################################################
# secondary_structure Tests 
###########################################################################
	add_executable(ss_pose_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/ss_pose_unittest.cpp)
	target_link_libraries(ss_pose_unittest secondary_structure_lib  )
	add_test(ss_pose_unittest ss_pose_unittest)
	
    add_executable(ss_chain_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/ss_chain_unittest.cpp)
	target_link_libraries(ss_chain_unittest secondary_structure_lib  )
	add_test(ss_chain_unittest ss_chain_unittest)
	
    add_executable(ss_structure_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/ss_structure_unittest.cpp)
	target_link_libraries(ss_structure_unittest secondary_structure_lib  )
	add_test(ss_structure_unittest ss_structure_unittest)
	
    add_executable(ss_residue_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/ss_residue_unittest.cpp)
	target_link_libraries(ss_residue_unittest secondary_structure_lib  )
	add_test(ss_residue_unittest ss_residue_unittest)
	
    add_executable(ss_motif_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/ss_motif_unittest.cpp)
	target_link_libraries(ss_motif_unittest secondary_structure_lib  )
	add_test(ss_motif_unittest ss_motif_unittest)
	
    add_executable(secondary_structure_parser_unittest ${RNAMAKE_UNITTEST_DIR}/secondary_structure/secondary_structure_parser_unittest.cpp)
	target_link_libraries(secondary_structure_parser_unittest secondary_structure_lib  )
	add_test(secondary_structure_parser_unittest secondary_structure_parser_unittest)
###########################################################################
# eternabot Tests 
###########################################################################
###########################################################################
# structure Tests 
###########################################################################
	add_executable(structure_unittest ${RNAMAKE_UNITTEST_DIR}/structure/structure_unittest.cpp)
	target_link_libraries(structure_unittest structure_lib  )
	add_test(structure_unittest structure_unittest)
	
    add_executable(basepair_state_unittest ${RNAMAKE_UNITTEST_DIR}/structure/basepair_state_unittest.cpp)
	target_link_libraries(basepair_state_unittest structure_lib  )
	add_test(basepair_state_unittest basepair_state_unittest)
	
    add_executable(pdb_parser_unittest ${RNAMAKE_UNITTEST_DIR}/structure/pdb_parser_unittest.cpp)
	target_link_libraries(pdb_parser_unittest structure_lib  )
	add_test(pdb_parser_unittest pdb_parser_unittest)
	
    add_executable(chain_closure_unittest ${RNAMAKE_UNITTEST_DIR}/structure/chain_closure_unittest.cpp)
	target_link_libraries(chain_closure_unittest structure_lib  )
	add_test(chain_closure_unittest chain_closure_unittest)
	
    add_executable(chain_unittest ${RNAMAKE_UNITTEST_DIR}/structure/chain_unittest.cpp)
	target_link_libraries(chain_unittest structure_lib  )
	add_test(chain_unittest chain_unittest)
	
    add_executable(atom_unittest ${RNAMAKE_UNITTEST_DIR}/structure/atom_unittest.cpp)
	target_link_libraries(atom_unittest structure_lib  )
	add_test(atom_unittest atom_unittest)
	
    add_executable(residue_type_unittest ${RNAMAKE_UNITTEST_DIR}/structure/residue_type_unittest.cpp)
	target_link_libraries(residue_type_unittest structure_lib  )
	add_test(residue_type_unittest residue_type_unittest)
	
    add_executable(basepair_unittest ${RNAMAKE_UNITTEST_DIR}/structure/basepair_unittest.cpp)
	target_link_libraries(basepair_unittest structure_lib  )
	add_test(basepair_unittest basepair_unittest)
	
    add_executable(residue_unittest ${RNAMAKE_UNITTEST_DIR}/structure/residue_unittest.cpp)
	target_link_libraries(residue_unittest structure_lib  )
	add_test(residue_unittest residue_unittest)
	
    add_executable(cif_parser_unittest ${RNAMAKE_UNITTEST_DIR}/structure/cif_parser_unittest.cpp)
	target_link_libraries(cif_parser_unittest structure_lib  )
	add_test(cif_parser_unittest cif_parser_unittest)
###########################################################################
# motif Tests 
###########################################################################
	add_executable(motif_state_ensemble_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_state_ensemble_unittest.cpp)
	target_link_libraries(motif_state_ensemble_unittest motif_lib  )
	add_test(motif_state_ensemble_unittest motif_state_ensemble_unittest)
	
    add_executable(motif_ensemble_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_ensemble_unittest.cpp)
	target_link_libraries(motif_ensemble_unittest motif_lib  )
	add_test(motif_ensemble_unittest motif_ensemble_unittest)
	
    add_executable(motif_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_unittest.cpp)
	target_link_libraries(motif_unittest motif_lib  )
	add_test(motif_unittest motif_unittest)
	
    add_executable(motif_factory_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_factory_unittest.cpp)
	target_link_libraries(motif_factory_unittest motif_lib  )
	add_test(motif_factory_unittest motif_factory_unittest)
	
    add_executable(motif_to_secondary_structure_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_to_secondary_structure_unittest.cpp)
	target_link_libraries(motif_to_secondary_structure_unittest motif_lib  )
	add_test(motif_to_secondary_structure_unittest motif_to_secondary_structure_unittest)
	
    add_executable(motif_scorer_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_scorer_unittest.cpp)
	target_link_libraries(motif_scorer_unittest motif_lib  )
	add_test(motif_scorer_unittest motif_scorer_unittest)
	
    add_executable(motif_state_unittest ${RNAMAKE_UNITTEST_DIR}/motif/motif_state_unittest.cpp)
	target_link_libraries(motif_state_unittest motif_lib  )
	add_test(motif_state_unittest motif_state_unittest)
###########################################################################
# motif_tools Tests 
###########################################################################
	add_executable(segmenter_unittest ${RNAMAKE_UNITTEST_DIR}/motif_tools/segmenter_unittest.cpp)
	target_link_libraries(segmenter_unittest motif_tools_lib  )
	add_test(segmenter_unittest segmenter_unittest)
###########################################################################
# resources Tests 
###########################################################################
	add_executable(added_motif_library_unittest ${RNAMAKE_UNITTEST_DIR}/resources/added_motif_library_unittest.cpp)
	target_link_libraries(added_motif_library_unittest resources_lib  )
	add_test(added_motif_library_unittest added_motif_library_unittest)
	
    add_executable(motif_sqlite_connection_unittest ${RNAMAKE_UNITTEST_DIR}/resources/motif_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_sqlite_connection_unittest resources_lib  )
	add_test(motif_sqlite_connection_unittest motif_sqlite_connection_unittest)
	
    add_executable(motif_sqlite_library_unittest ${RNAMAKE_UNITTEST_DIR}/resources/motif_sqlite_library_unittest.cpp)
	target_link_libraries(motif_sqlite_library_unittest resources_lib  )
	add_test(motif_sqlite_library_unittest motif_sqlite_library_unittest)
	
    add_executable(motif_state_sqlite_library_unittest ${RNAMAKE_UNITTEST_DIR}/resources/motif_state_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_sqlite_library_unittest resources_lib  )
	add_test(motif_state_sqlite_library_unittest motif_state_sqlite_library_unittest)
	
    add_executable(resource_mananger_unittest ${RNAMAKE_UNITTEST_DIR}/resources/resource_mananger_unittest.cpp)
	target_link_libraries(resource_mananger_unittest resources_lib  )
	add_test(resource_mananger_unittest resource_mananger_unittest)
	
    add_executable(motif_ensemble_sqlite_connection_unittest ${RNAMAKE_UNITTEST_DIR}/resources/motif_ensemble_sqlite_connection_unittest.cpp)
	target_link_libraries(motif_ensemble_sqlite_connection_unittest resources_lib  )
	add_test(motif_ensemble_sqlite_connection_unittest motif_ensemble_sqlite_connection_unittest)
	
    add_executable(motif_state_ensemble_sqlite_library_unittest ${RNAMAKE_UNITTEST_DIR}/resources/motif_state_ensemble_sqlite_library_unittest.cpp)
	target_link_libraries(motif_state_ensemble_sqlite_library_unittest resources_lib  )
	add_test(motif_state_ensemble_sqlite_library_unittest motif_state_ensemble_sqlite_library_unittest)
###########################################################################
# motif_data_structure Tests 
###########################################################################
	add_executable(motif_state_graph_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_state_graph_unittest.cpp)
	target_link_libraries(motif_state_graph_unittest motif_data_structure_lib  )
	add_test(motif_state_graph_unittest motif_state_graph_unittest)
	
    add_executable(motif_state_tree_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_state_tree_unittest.cpp)
	target_link_libraries(motif_state_tree_unittest motif_data_structure_lib  )
	add_test(motif_state_tree_unittest motif_state_tree_unittest)
	
    add_executable(motif_topology_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_topology_unittest.cpp)
	target_link_libraries(motif_topology_unittest motif_data_structure_lib  )
	add_test(motif_topology_unittest motif_topology_unittest)
	
    add_executable(motif_tree_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_tree_unittest.cpp)
	target_link_libraries(motif_tree_unittest motif_data_structure_lib  )
	add_test(motif_tree_unittest motif_tree_unittest)
	
    add_executable(motif_state_ensemble_tree_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_state_ensemble_tree_unittest.cpp)
	target_link_libraries(motif_state_ensemble_tree_unittest motif_data_structure_lib  )
	add_test(motif_state_ensemble_tree_unittest motif_state_ensemble_tree_unittest)
	
    add_executable(motif_merger_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_merger_unittest.cpp)
	target_link_libraries(motif_merger_unittest motif_data_structure_lib  )
	add_test(motif_merger_unittest motif_merger_unittest)
	
    add_executable(motif_graph_unittest ${RNAMAKE_UNITTEST_DIR}/motif_data_structure/motif_graph_unittest.cpp)
	target_link_libraries(motif_graph_unittest motif_data_structure_lib  )
	add_test(motif_graph_unittest motif_graph_unittest)
###########################################################################
# thermo_fluctuation Tests 
###########################################################################
	add_executable(thermo_fluc_sampler_unittest ${RNAMAKE_UNITTEST_DIR}/thermo_fluctuation/thermo_fluc_sampler_unittest.cpp)
	target_link_libraries(thermo_fluc_sampler_unittest thermo_fluctuation_lib  )
	add_test(thermo_fluc_sampler_unittest thermo_fluc_sampler_unittest)
	
    add_executable(thermo_fluc_simulation_unittest ${RNAMAKE_UNITTEST_DIR}/thermo_fluctuation/thermo_fluc_simulation_unittest.cpp)
	target_link_libraries(thermo_fluc_simulation_unittest thermo_fluctuation_lib  )
	add_test(thermo_fluc_simulation_unittest thermo_fluc_simulation_unittest)
###########################################################################
# motif_search Tests 
###########################################################################
	add_executable(exhaustive_search_unittest ${RNAMAKE_UNITTEST_DIR}/motif_search/exhaustive_search_unittest.cpp)
	target_link_libraries(exhaustive_search_unittest motif_search_lib  )
	add_test(exhaustive_search_unittest exhaustive_search_unittest)
	
    add_executable(monte_carlo_search_unittest ${RNAMAKE_UNITTEST_DIR}/motif_search/monte_carlo_search_unittest.cpp)
	target_link_libraries(monte_carlo_search_unittest motif_search_lib  )
	add_test(monte_carlo_search_unittest monte_carlo_search_unittest)
	
    add_executable(motif_state_search_unittest ${RNAMAKE_UNITTEST_DIR}/motif_search/motif_state_search_unittest.cpp)
	target_link_libraries(motif_state_search_unittest motif_search_lib  )
	add_test(motif_state_search_unittest motif_state_search_unittest)
###########################################################################
# sequence_optimization Tests 
###########################################################################
	add_executable(sequence_optimizer_unittest ${RNAMAKE_UNITTEST_DIR}/sequence_optimization/sequence_optimizer_unittest.cpp)
	target_link_libraries(sequence_optimizer_unittest sequence_optimization_lib  )
	add_test(sequence_optimizer_unittest sequence_optimizer_unittest)
####################################################################################################

