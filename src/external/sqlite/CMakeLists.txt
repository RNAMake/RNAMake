set( SQLITE_DIR ${RNAMAKE_EXT_SRC_DIR}/sqlite/ )
include_directories( ${SQLITE_DIR} )

set( sql_src ${SQLITE_DIR}/sqlite3.c )

set_source_files_properties( ${SQLITE_DIR}/sqlite3.h
            PROPERTIES
            LANGUAGE CXX
            )

set_source_files_properties(
            sql_src
            PROPERTIES
            LANGUAGE C
            LINKER CXX_STATIC_LIBRARY_LINKER
            )


if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_library(sqlite3 STATIC ${SQLITE_DIR}/sqlite3.c )
else()
    add_library(sqlite3 ${SQLITE_DIR}/sqlite3.c )
endif()

target_link_libraries(sqlite3 -ldl )
