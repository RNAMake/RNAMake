include_directories(${RNAMAKE}/src/external/sqlite/)

set(sql_src
        ${RNAMAKE}/src/external/sqlite/sqlite3.c
        )

set_source_files_properties(
        ${RNAMAKE}/src/external/sqlite/sqlite3.h
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
    add_library(sqlite3 STATIC ${RNAMAKE}/src/external/sqlite/sqlite3.c )
    target_link_libraries(sqlite3 -static -ldl)

else()

    add_library(sqlite3 ${RNAMAKE}/src/external/sqlite/sqlite3.c )
    target_link_libraries(sqlite3)



endif()
