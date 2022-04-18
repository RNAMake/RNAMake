//
//  settings.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <fstream>

// RNAMake Headers
#include <base/exception.hpp>
#include <base/types.hpp>

namespace base::path {

String filename(String const &str_path) {
  Path path = {str_path};
  if (!path.has_filename()) {
    String msg = str_path + " does not have a filename!";
    base::log_and_throw<base::InputException>(msg);
  }
  return path.filename();
}

String parent_dir(String const &str_path) {
  Path path = {str_path};
  if (!path.has_parent_path()) {
    String msg = str_path + " does not have a parent path!";
    base::log_and_throw<base::InputException>(msg);
  }
  return path.parent_path();
}

String rnamake_path() {
  char *base_path = std::getenv("RNAMAKE");
  if (base_path == nullptr) {
    String msg = "cannot find environmental path RNAMAKE, please set it"
                 "should be set to $PATH/RNAMake where $PATH is where you "
                 "installed RNAMake";
    base::log_and_throw<base::ResourceException>(msg);
  }
  if (!std::filesystem::exists(base_path)) {
    String msg = "environmental path RNAMAKE is set but is not a valid path";
    base::log_and_throw<base::ResourceException>(msg);
  }
  String path = {base_path};
  if (!std::filesystem::exists(path + "/resources")) {
    String msg = "environmental path RNAMAKE is set but does not contain a "
                 "resources directory!";
    base::log_and_throw<base::ResourceException>(msg);
  }
  return path;
}

String resources_path() { return rnamake_path() + "resources/"; }

String unittest_resource_path() {
  return rnamake_path() + "unittests/unittest_resources/";
}

void get_lines_from_file(String const & fname, Strings & lines) {
  if(!std::filesystem::exists(fname)) {
    String msg = "The file " + fname + " does not exist!";
    base::log_and_throw<base::InputException>(msg);
  }

  String line;
  std::ifstream input;
  input.open(fname);
  while (input.good()) {
    getline(input, line);
    lines.push_back(line);
  }
  lines.pop_back();
  if(lines.empty()) {
    LOG_WARNING << "there are no lines in " << fname;
  }
  input.close();
}

} // namespace base::path
