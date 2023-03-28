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
  auto seg2 = rm.get_segment(resource_management::SegmentInfo{"HELIX.IDEAL.4"});
  sg.add(seg1);
  sg.add(seg2, 0, sg.get_end_name(0, 1));
  auto persistence = Persistence();
  ofstream file;
  file.open("test_dir/canary_file");
  std::cout << "Running test cases\n";
  persistence.save_to_database(sg, "test", "test_seg");

  // Actual test cases
  SUBCASE("General testing") {
    bool db_dir = filesystem::is_directory("test_dir");
    CHECK(db_dir == true);
  }

  SUBCASE("Does not override existing directory") {
    bool canary_still_exists = filesystem::exists("test_dir/canary_file");
    CHECK(canary_still_exists == true);
  }

  // SUBCASE("Creates new database") {
  //   // Database was created in previous tests
  //   bool db_exists = filesystem::exists("user_database_dir/user_database.db");
  //   CHECK(db_exists == true);
  // }

  // SUBCASE("Writes data to database") {
  //   SUBCASE("Does not write existing segment to databse") {
  //     // Write segments of sg to database
  //     // Ensure the segment table has 2 records
  //     // Create an identical segment graph object called sg2
  //     // Write segments of sg2 to databse
  //     // Ensure the segment table has only 2 records
  //   }

  //   SUBCASE("Writes new segment and doesn't write old segment") {
  //     // Write sg to databse
  //     // Ensure the segment table has 2 records
  //     // Create sg2 object with seg1 and new segment
  //     // Write sg2 to database
  //     // Ensure the segment table has only 3 records
  //   }
  // }

  // SUBCASE("Does not override existing database") {

  // }
}
