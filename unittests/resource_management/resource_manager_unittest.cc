//
// Created by Joseph Yesselman on 1/3/18.
//

#include "../common.hpp"
#include <iostream>

#include <resource_management/resource_manager.h>

using namespace resource_management;

TEST_CASE("Test resource manager") {
  // init_unittest_safe_logging();

  ResourceManager rm;
  // okay !
  auto const &rm_ref = rm;
  // not okay
  // auto rm1 = rm;

  // get segment
  auto seg = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
  CHECK(seg.get_name() == "HELIX.IDEAL.2");
}