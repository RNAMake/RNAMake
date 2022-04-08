//
//  settings.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

// RNAMake Headers
#include <base/exception.hpp>
#include <base/types.hpp>

namespace base::path {

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

/*String resources_path() {
  String base_path = rnamake_path();
  return base_path + "/resources/";
}

String unittest_resource_path() {
  return rnamake_path() + "/unittests/unittest_resources/";
}      */

} // namespace base::path
