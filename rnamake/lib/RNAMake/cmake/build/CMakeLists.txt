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

add_executable(unit ../../unittests/main.cc)
target_link_libraries(unit ${link_libraries})
#set_target_properties(unit PROPERTIES COMPILE_FLAGS “${COMPILE_FLAGS}”)
add_custom_target(unit_symlink ALL)
add_custom_command(TARGET unit_symlink POST_BUILD COMMAND cmake -E create_symlink ../../bin/unittests/unit unit)