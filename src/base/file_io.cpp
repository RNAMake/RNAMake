//
//  FileIO.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/29/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//


#include "base/file_io.hpp"

namespace base {

void get_lines_from_file(String const fname, Strings & lines) {
  if (!file_exists(fname)) {
    throw("ERROR: The file " + fname + " does not exist");
  }
  String line;
  std::ifstream input;
  input.open(fname);
  while (input.good()) {
    getline(input, line);
    lines.push_back(line);
  }
  input.close();
}

}  // namespace base