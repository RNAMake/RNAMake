//
// Created by Joseph Yesselman on 1/3/18.
//

#include "../common.hpp"
#include <iostream>

#include <resource_management/resource_manager.h>

TEST_CASE("Test resource manager") {
  // init_unittest_safe_logging();

  resource_management::ResourceManager rm;
  // okay !
  auto const &rm_ref = rm;
  // not okay
  // auto rm1 = rm;

  // get segment
  // REQUIRE_NOTHROW(rm.get_segment(StringStringMap{{"name", "HELIX.IDEAL"
  //                                                          ".2"}}));
}