//
// Created by Erik Whiting on 2/13/23
//

#include "../common.hpp"
#include <iostream>

#include <structure/base/pose.hpp>
#include <structure/base/segment.hpp>
#include <resource_management/resource_manager.h>

using namespace structure::base;
using namespace resource_management;

TEST_CASE("Test Segment") {
  SUBCASE("Equality operator") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    auto seg2 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    CHECK(*seg1 == *seg2);
  }

  SUBCASE("Inequality operator") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    auto seg2 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.4"});
    CHECK(*seg1 != *seg2);
  }
}
