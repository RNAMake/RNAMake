
//
// Created by Chris Jurich  on 4/16/20.
//

#ifndef BASE_LIB_BACKTRACE_H
#define BASE_LIB_BACKTRACE_H


#include <iostream>
#include <stdio.h>
#include <execinfo.h>
#include <cxxabi.h>


std::string
demangle( std::string );

void
print_backtrace();

#endif //BASE_LIB_BACKTRACE_H
