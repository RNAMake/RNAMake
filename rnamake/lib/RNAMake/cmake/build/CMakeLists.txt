cmake_minimum_required(VERSION 2.8)

project (test)
set(CMAKE_BUILD_TYPE Release)

include(compiler.cmake)

include_directories(../../src/)
include_directories(../../unittests/)

include(libraries.cmake )

set( previous_library )
set( link_libraries ${previous_library} )
find_library(SQLITE3_LIBRARY sqlite3)
foreach( library ${libraries} )
	include(${library}.cmake)
	add_library(${library} SHARED ${${library}_files})
    target_link_libraries( ${library} ${previous_library} ${SQLITE3_LIBRARY} )
	set( link_libraries ${library};${link_libraries} )
	set( previous_library ${library} )
endforeach()

include(unittest_apps.cmake)
include(apps.cmake )

