set(RNAMAKE_INTEGRATION_DIR ${RNAMAKE}/integration_tests/ )


enable_testing()

include_directories( ${RNAMAKE}/src/ )
include_directories( ${RNAMAKE}/apps/ )

include_directories( ${RNAMAKE_INTEGRATION_DIR} )

if ( NOT EXISTS ${BUILD_DIR}/integration_tests )
    make_directory( ${BUILD_DIR}/integration_tests )
endif()

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BUILD_DIR}/integration_tests )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BUILD_DIR}/integration_tests )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BUILD_DIR}/integration_tests )


#####################################################################################################
## Integration tests Declarations 
#####################################################################################################
	add_executable(TTR_integration ${RNAMAKE_INTEGRATION_DIR}/TTR/TTR_integration.cpp)
	target_link_libraries(TTR_integration RNAMake  )
	add_test(TTR_integration TTR_integration)
####################################################################################################
	add_executable(rebuild_p4_p6_integration ${RNAMAKE_INTEGRATION_DIR}/REBUILD_P4_P6/rebuild_p4_p6_integration.cpp)
	target_link_libraries(rebuild_p4_p6_integration RNAMake  )
	add_test(rebuild_p4_p6_integration rebuild_p4_p6_integration)
####################################################################################################
	add_executable(build_around_protein ${RNAMAKE_INTEGRATION_DIR}/BUILD_AROUND_PROTEIN/build_around_protein.cpp)
	target_link_libraries(build_around_protein RNAMake  )
	add_test(build_around_protein build_around_protein)
####################################################################################################
#    set_source_files_properties(  ${RNAMAKE_INTEGRATION_DIR}/tools.h 
#        PROPERTIES
#        LANGUAGE CXX
#        )
#	add_executable(tools ${RNAMAKE_INTEGRATION_DIR}/tools.h)
#	target_link_libraries(tools RNAMake  )
#	add_test(tools tools)
####################################################################################################
	add_executable(TTR_supplied_motif_integration ${RNAMAKE_INTEGRATION_DIR}/TTR_SUPPLIED_MOTIF/TTR_supplied_motif_integration.cpp)
	target_link_libraries(TTR_supplied_motif_integration RNAMake  )
	add_test(TTR_supplied_motif_integration TTR_supplied_motif_integration)
####################################################################################################
	add_executable(TTR_MC_integration ${RNAMAKE_INTEGRATION_DIR}/TTR_MC/TTR_MC_integration.cpp)
	target_link_libraries(TTR_MC_integration RNAMake  )
	add_test(TTR_MC_integration TTR_MC_integration)
####################################################################################################

