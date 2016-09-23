cmake_minimum_required(VERSION 2.8)

project (test)
set(CMAKE_BUILD_TYPE Release)

include(compiler.cmake)

include_directories(../../src/)
include_directories(../../unittests/)
include_directories(../../unittests_new/)

find_library(SQLITE3_LIBRARY sqlite3)

include(build.cmake)
include(apps.cmake)
