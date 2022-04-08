//
//  FileIO.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__FileIO__
#define __REDESIGNC__FileIO__

///////////////////////////////////////////////////////////////////////////////
// standard includes                                                         //
///////////////////////////////////////////////////////////////////////////////

#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
// rnamake include s                                                         //
///////////////////////////////////////////////////////////////////////////////

#include <base/types.hpp>

// This all needs to be replaced by https://en.cppreference.com/w/cpp/filesystem/path
namespace base {


void get_lines_from_file(String, Strings & /* return */) noexcept(false);

}  // namespace base

#endif /* defined(__REDESIGNC__FileIO__) */
