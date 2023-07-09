
//
// Created by Chris Jurich  on 4/16/20.
//

#ifndef BASE_LIB_BACKTRACE_H
#define BASE_LIB_BACKTRACE_H

#include <base/file_io.h>
#include <base/log.h>
#include <base/types.h>
#include <stdio.h>

#include <fstream>
#include <iostream>

#if defined(_WIN32) || defined(_WIN64)

#include <dbghelp.h>
#include <windows.h>

#else
#include <execinfo.h>
#endif

#include <cxxabi.h>

namespace base {

std::string demangle(std::string);

void print_backtrace();

void save_backtrace();

}  // namespace base
#endif  // BASE_LIB_BACKTRACE_H
