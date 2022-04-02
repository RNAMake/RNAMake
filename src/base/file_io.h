//
//  FileIO.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__FileIO__
#define __REDESIGNC__FileIO__

#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// RNAMake Headers
#include <base/exception.h>

#include "base/types.h"

namespace base {

inline bool file_exists(String const& name) {
  auto testfile = std::ifstream(name);
  const auto good = testfile.is_open();
  testfile.close();
  return good;
  // struct stat buffer;
  // return (stat(name.c_str(), &buffer) == 0);
}

inline int is_dir(String const& path) {
  struct stat info;

  if (stat(path.c_str(), &info) != 0) {
    return 0;
  } else if (info.st_mode & S_IFDIR) {
    return 1;
  } else {
    return 0;
  }

  //    printf( "cannot access %s\n", pathname );
  // else if( info.st_mode & S_IFDIR )  // S_ISDIR() doesn't exist on my windows
  //    printf( "%s is a directory\n", pathname );
  // else
  // printf( "%s is no directory\n", pathname );
  //
  // if (!file_exists(path)) { return 0; }
  // struct stat st;
  // lstat(path.c_str(), &st);
  // if (S_ISDIR(st.st_mode)) {
  //    return 1;
  //} else {
  //    return 0;
  //}
}

Strings get_lines_from_file(String) noexcept(false);

}  // namespace base

#endif /* defined(__REDESIGNC__FileIO__) */
