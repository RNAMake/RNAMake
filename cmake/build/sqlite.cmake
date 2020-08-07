include_directories(../../third-party/sqlite/)

set(sql_src
        ../../third-party/sqlite/sqlite3.c
        )

set_source_files_properties(
            ../../third-party/sqlite.h
            PROPERTIES
            LANGUAGE CXX
            )

set_source_files_properties(
            sql_src
            PROPERTIES
            LANGUAGE C
            LINKER CXX_STATIC_LIBRARY_LINKER
            )


add_library(sqlite3 ../../third-party/sqlite/sqlite3.c )


target_link_libraries(sqlite3)
