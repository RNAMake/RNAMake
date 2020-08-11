include_directories(../../src/external/sqlite/)

set(sql_src
        ../../src/external/sqlite/sqlite3.c
        )

set_source_files_properties(
        ../../src/external/sqlite/sqlite3.h
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
    add_library(sqlite3 STATIC ../../src/external/sqlite/sqlite3.c )
    target_link_libraries(sqlite3)

else()
add_library(sqlite3 ../../src/external/sqlite/sqlite3.c )
target_link_libraries(sqlite3)



endif()
