
//
// Created by Joseph Yesselman on 11/12/17.
//

#ifndef RNAMAKE_NEW_PATHS_H
#define RNAMAKE_NEW_PATHS_H

#include <base/types.h>

namespace base {

  String
  get_os_name();

  String
  base_path();

  String
  resources_path();

  String
  unittest_resources_path();

}

#endif //RNAMAKE_NEW_PATHS_H