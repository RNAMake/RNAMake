include_directories(../../src/external/nlohmann/)

set(json_src
        ../../src/external/nlohmann/json.hpp
        )

    set_source_files_properties(
            ../../src/external/nlohmann/json.hpp
                PROPERTIES
                LANGUAGE CXX
                LINKER_LANGUAGE CXX
                )
    
    
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    #add_library(sqlite3 STATIC ../../src/external/sqlite/sqlite3.c )
    #target_link_libraries(sqlite3 -static)
    
else()

    add_library(json ${json_src})
    set_target_properties(json 
                            PROPERTIES
                            LANGUAGE CXX
                            )
    target_precompile_headers(json 
            PUBLIC
            ../../src/external/nlohmann/json.hpp
        ) 
    target_link_libraries(json)



endif()
