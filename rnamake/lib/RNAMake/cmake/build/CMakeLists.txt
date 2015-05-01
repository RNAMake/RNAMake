cmake_minimum_required(VERSION 2.8)
CMAKE_POLICY( SET CMP0042 NEW )

project (test)
set(CMAKE_BUILD_TYPE Release)

include(compiler.cmake)

include_directories(../../src/)
include_directories(../../unittests/)

include(libraries.cmake )

set( previous_library )
set( link_libraries ${previous_library} )
foreach( library ${libraries} )
	include(${library}.cmake)
	add_library(${library} SHARED ${${library}_files})
	target_link_libraries( ${library} ${previous_library} )
	set( link_libraries ${library};${link_libraries} )
	set( previous_library ${library} )
endforeach()

find_library(SQLITE3_LIBRARY sqlite3)
add_executable(all_tests ../../unittests/all_tests.cc)
target_link_libraries(all_tests ${link_libraries} ${SQLITE3_LIBRARY})
#set_target_properties(all_tests PROPERTIES COMPILE_FLAGS “${COMPILE_FLAGS}”)
add_custom_target(unit_symlink ALL)
add_custom_command(TARGET unit_symlink POST_BUILD COMMAND cmake -E create_symlink ../../bin/unittests/all_tests all_tests)