//
//  settings.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>

// RNAMake Headers
#include <base/log.hpp>
#include <base/types.hpp>

namespace base {

String get_os_name() {
#if defined(_WIN32) || defined(_WIN64)
  return String("Windows");
#elif defined(__unix) || defined(__unix__)
  return String("linux");
#elif defined(__APPLE__)
  return String("OSX");
#elif defined(__linux__)
  return String("linux");
#else
#warning \
    "Could not deduce operating system type. You WILL encounter a runtime error."

  throw std::runtime_error("cannot determine operating system");
#endif
}

String base_dir() {
  char *base_path;
  base_path = std::getenv("RNAMAKE");
  if (base_path == NULL) {
    LOG_ERROR << "cannot find environmental path RNAMAKE, please set it";
    LOG_ERROR << "should be set to $PATH/RNAMake where $PATH is where you "
                 "installed RNAMake";
    exit(EXIT_FAILURE);
  }
  return String(base_path);
}

String resources_path() {
  String base_path = base_dir();
  return base_path + "/resources/";
}

String lib_path() {
  String base_path = base_dir();
  return base_path;
}

String motif_dirs() {
  String base_path = base_dir();
  return base_path + "/resources/motifs/";
}

String unittest_resource_dir() {
  return base_dir() + "/unittests/unittest_resources/";
}

}  // namespace base
