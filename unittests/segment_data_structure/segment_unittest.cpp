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
  SUBCASE("Secondary Structure to_str Method") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    // Wait, what?
    String s = seg1->to_str();
    std::cout << s << std::endl;
  }
}
