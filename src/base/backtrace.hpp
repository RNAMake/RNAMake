
//
// Created by Chris Jurich  on 4/16/20.
//

#ifndef BASE_LIB_BACKTRACE_H
#define BASE_LIB_BACKTRACE_H

///////////////////////////////////////////////////////////////////////////////
// standard includes                                                         //
///////////////////////////////////////////////////////////////////////////////

#if defined(_WIN32) || defined(_WIN64)

#include <dbghelp.h>
#include <windows.h>

#else
#include <execinfo.h>
#endif

#include <cxxabi.h>
#include <cstdio>

#include <fstream>
#include <iostream>

// rnamake includes
#include <base/log.hpp>
#include <base/types.hpp>

namespace base {

String demangle(String);

void print_backtrace();

void save_backtrace();

}  // namespace base
#endif  // BASE_LIB_BACKTRACE_H
