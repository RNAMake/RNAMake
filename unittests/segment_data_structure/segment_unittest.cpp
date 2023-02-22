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
  SUBCASE("Segmet to_str Method") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    String s = seg1->to_str();
    std::cout << s << "\n\n";
  }

  SUBCASE("Secondary structure string method") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    String s = seg1->secondary_structure_to_str();
    std::cout << s << std::endl;
  }

  SUBCASE("to_str is true to get_str") {
    ResourceManager rm;
    auto seg1 = rm.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    String segment_string = seg1->to_str();
    auto seg2 = structure::all_atom::get_segment_from_str(segment_string);
    // May need to implement an `operator==` method
    // CHECK(seg1 == seg2);
  }
}
