
set(RNAMAKE_APP_DIR ${RNAMAKE}/apps/ )

include_directories( ${RNAMAKE_APP_DIR} )

if ( NOT EXISTS ${BUILD_DIR}/exe )
    make_directory( ${BUILD_DIR}/exe )
endif()

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BUILD_DIR}/exe )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD_DIR}/exe )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}/exe )

###########################################################################
# ddg_tecto
###########################################################################
	add_executable( ddg_tecto ${RNAMAKE_APP_DIR}/ddg_tecto/ddg_tecto.cc)
	target_link_libraries(ddg_tecto RNAMake   )
###########################################################################
# build_sqlite_libraries
###########################################################################
	add_executable( build_sqlite_libraries ${RNAMAKE_APP_DIR}/build_sqlite_libraries/build_sqlite_libraries.cpp)
	target_link_libraries(build_sqlite_libraries RNAMake   )
###########################################################################
# design_rna_scaffold
###########################################################################
add_executable( design_rna_scaffold ${RNAMAKE_APP_DIR}/design_rna_scaffold/design_rna_scaffold.cpp ${RNAMAKE_APP_DIR}/design_rna_scaffold/main.cpp)
	target_link_libraries(design_rna_scaffold RNAMake   )
###########################################################################
# build_motif_graph
###########################################################################
	add_executable( build_motif_graph ${RNAMAKE_APP_DIR}/build_motif_graph/build_motif_graph.cxx)
	target_link_libraries(build_motif_graph RNAMake   )
###########################################################################
# compare_db_files
###########################################################################
	add_executable( compare_db_files ${RNAMAKE_APP_DIR}/compare_db_files/compare_db_files.cpp)
	target_link_libraries(compare_db_files RNAMake   )
###########################################################################
# sample_helix
###########################################################################
	add_executable( sample_helix ${RNAMAKE_APP_DIR}/sample_helix/sample_helix.cpp)
	target_link_libraries(sample_helix RNAMake   )
####################################################################################################
