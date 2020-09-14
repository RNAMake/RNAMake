
//
// Created by Chris Jurich  on 4/16/20.
//

#ifndef BASE_LIB_BACKTRACE_H
#define BASE_LIB_BACKTRACE_H

#include <fstream>
#include <iostream>
#include <stdio.h>

#include <base/log.h>
#include <base/types.h>
#include <base/file_io.h>

#if defined(_WIN32) || defined(_WIN64)
#   warning RNAMake backtrace is not currently supported on windows platforms
#else
#   include <execinfo.h>
#   include <cxxabi.h>
#endif



namespace base {

std::string
demangle( std::string );

void
print_backtrace();

void
save_backtrace();

}
#endif //BASE_LIB_BACKTRACE_H
