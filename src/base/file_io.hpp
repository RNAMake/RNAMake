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

namespace base {

inline bool file_exists(String const& name) {
  auto testfile = std::ifstream(name);
  bool good = testfile.is_open();
  testfile.close();
  return good;
}

inline int is_dir(String const& path) {
  struct stat info{};
  if (stat(path.c_str(), &info) != 0) {
    return 0;
  } else if (info.st_mode & S_IFDIR) {
    return 1;
  } else {
    return 0;
  }
}

void get_lines_from_file(String, Strings & /* return */) noexcept(false);

}  // namespace base

#endif /* defined(__REDESIGNC__FileIO__) */
