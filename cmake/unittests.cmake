###########################################################################
# base Tests
###########################################################################
foreach(component
        cl_option;
        env_manager;
        application;
        log;
        string;
        option
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/base/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest base_lib)
endforeach()
###########################################################################
# math Tests
###########################################################################
foreach(component
        rotation;
        xyz_matrix;
        hashing;
        quaternion)
    add_executable(${component}_unittest ${RNAMAKE}/unittests/math/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest math_lib)
endforeach()
###########################################################################
# data_structure Tests
###########################################################################
foreach(component
        graph;
        graph_iter_list;
        tree;
        graph_adjacency_list
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/data_structure/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest data_structure_lib)
endforeach()
###########################################################################
# util Tests
###########################################################################
foreach(component
        x3dna;
        uuid;
        sqlite3_connection;
        csv;
        steric_lookup;
        dssr
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/util/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest util_lib)
endforeach()
###########################################################################
# vienna Tests
###########################################################################
###########################################################################
# secondary_structure Tests
###########################################################################
foreach(component
        ss_chain;
        ss_motif;
        ss_structure;
        secondary_structure_parser;
        ss_pose;
        ss_residue
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/secondary_structure/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest secondary_structure_lib)
endforeach()
###########################################################################
# eternabot Tests
###########################################################################
###########################################################################
# structure Tests
###########################################################################
foreach(component
        chain;
        structure;
        basepair;
        residue_type;
        chain_closure;
        cif_parser;
        atom;
        basepair_state;
        pdb_parser;
        residue;
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/structure/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest structure_lib  )
endforeach()
###########################################################################
# motif Tests
###########################################################################
foreach(component
        motif_state_ensemble;
        motif_state;
        motif;
        motif_factory;
        motif_scorer;
        motif_to_secondary_structure
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/motif/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest motif_lib)
endforeach()
###########################################################################
# motif_tools Tests
###########################################################################
foreach(component
        segmenter
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/motif_tools/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest motif_tools_lib )
endforeach()
###########################################################################
# resources Tests
###########################################################################
foreach(component	
        resource_mananger;
	    motif_sqlite_connection;
	    motif_ensemble_sqlite_connection;
	    motif_state_sqlite_library;
	    added_motif_library;
	    motif_sqlite_library
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/resources/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest resources_lib  )
endforeach()
###########################################################################
# motif_data_structure Tests
###########################################################################
foreach(component
	    motif_topology;
	    motif_state_graph;
	    motif_merger;
	    motif_tree;
	    motif_state_tree;
	    motif_graph
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/motif_data_structure/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest motif_data_structure_lib  )

endforeach()
###########################################################################
# thermo_fluctuation Tests
###########################################################################
foreach(component	
        thermo_fluc_simulation;
        thermo_fluc_sampler
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/thermo_fluctuation/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest thermo_fluctuation_lib )
endforeach()
###########################################################################
# motif_search Tests
###########################################################################
foreach(component	
        motif_state_search;
	    monte_carlo_search;
	    exhaustive_search
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/motif_search/${component}_unittest.cpp)
    target_link_libraries(${component}_unittest motif_search_lib)
endforeach()
###########################################################################
# sequence_optimization Tests
###########################################################################
foreach(component
        sequence_optimizer
        )
    add_executable(${component}_unittest ${RNAMAKE}/unittests/sequence_optimization/${component}_unittest.cpp)
	target_link_libraries(${component}_unittest sequence_optimization_lib  )
endforeach()
