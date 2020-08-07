
//
// Created by Chris Jurich  on 4/16/20.
//

#ifndef BASE_LIB_BACKTRACE_H
#define BASE_LIB_BACKTRACE_H


#include <iostream>
#include <stdio.h>

#if defined(_WIN32) || defined(_WIN64)

#include <windows.h>
#include <DbgHelp.h>

#else
#include <execinfo.h>
#endif

#include <cxxabi.h>
namespace base {

std::string
demangle( std::string );

void
print_backtrace();
}
#endif //BASE_LIB_BACKTRACE_H
