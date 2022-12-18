//
// Created by Erik Whiting on 12/10/22.
//

#include "../common.hpp"
#include <iostream>
#include <filesystem>
#include <sqlite3.h>

#include "resource_management/resource_manager.h"
#include "segment_data_structure/segment_graph.h"
#include "segment_data_structure/persistence.h"

using namespace persistence;
using namespace std;

TEST_CASE("Test graph persistence") {
  // Test setup stuff
  resource_management::ResourceManager rm;
  segment_data_structure::SegmentGraphAllAtom sg;
  auto seg1 = rm.get_segment(resource_management::SegmentInfo{"HELIX.IDEAL.2"});
  auto seg2 = rm.get_segment(resource_management::SegmentInfo{"HELIX.IDEAL.2"});
  sg.add(seg1);
  sg.add(seg2, 0, sg.get_end_name(0, 1));

  // Actual test cases
  SUBCASE("Creates new directory") {
    Persistence::save_to_database(sg);
    bool db_dir = filesystem::is_directory("user_database_dir");
    CHECK(db_dir == true);
  }

  SUBCASE("Does not override existing directory") {
    ofstream file;
    file.open("user_database_dir/canary_file");
    Persistence::save_to_database(sg);
    bool canary_still_exists = filesystem::exists("user_database_dir/canary_file");
    CHECK(canary_still_exists == true);
  }

  SUBCASE("Creates new database") {
    // Database was created in previous tests
    bool db_exists = filesystem::exists("user_database_dir/user_database.db");
    CHECK(db_exists == true);
  }

  SUBCASE("Writes data to database") {

  }

  SUBCASE("Does not override existing database") {

  }
}
