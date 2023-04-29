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
  // When testing the object out of the user database, ignore
  // the difference in center when testing the segment objects
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
  std::cout << "Running test cases\n\n\n";

  // Actual test cases
  SUBCASE("General testing") {
    persistence.save_to_database(sg, "test", "test_seg");
    bool db_dir = filesystem::is_directory("test_dir");
    CHECK(db_dir == true);
  }

  SUBCASE("Does not override existing directory") {
    bool canary_still_exists = filesystem::exists("test_dir/canary_file");
    CHECK(canary_still_exists == true);
  }

  SUBCASE("Builds segment object from database") {
    std::cout << "About to run retireve_segment_from_database method\n";
    CHECK_NOTHROW(
      persistence.retrieve_segment_from_database("HELIX.IDEAL.2", "test_dir/test.db")
    );
    auto test_seg = persistence.retrieve_segment_from_database("HELIX.IDEAL.2", "test_dir/test.db");
  }

  SUBCASE("Builds segment graph object from databse") {
    std::cout << "About to run retrive_segment_graph_from_database method\n";
    CHECK_NOTHROW(
      persistence.retrieve_segment_graph_from_database("test_seg", "test_dir/test.db")
    );
  }

  SUBCASE("Check db/memory fidelity") {
    SUBCASE("Segment fidelity") {
      auto db_helix_ideal_2 = persistence.retrieve_segment_from_database("HELIX.IDEAL.2", "test_dir/test.db");
      auto mm_helix_ideal_2 = rm.get_segment(resource_management::SegmentInfo{"HELIX.IDEAL.2"});
      CHECK(*db_helix_ideal_2 == *mm_helix_ideal_2);
    }

    SUBCASE("Segment graph fidelity") {
      const segment_data_structure::SegmentGraphAllAtom db_graph = persistence.retrieve_segment_graph_from_database("test_seg", "test_dir/test.db");
      CHECK(sg == db_graph);
    }
  }

  // Tests!!!!
  // Obvious: build like ten segments in a row, and build the exact same way
  // This means checking that the coordinates and bp ref frames are the same
  // get_ref_frame, it's a 3x3 matrix (fairly large number so aggregate error isn't bad),
  // then for each motif object, save it to databse, rebuild it as copy,
  // compare original to the one saved. For each segment, check each of its ends
  // so the origin (make function to check to segment graphs).
  // check if bp ref frames are very close to each other and that origins (center) are simlar
  // Check atomic coordinates in segment. Do this by ...
  // Gonna have to think of a way to
  // Each residue has atoms, want atom coordinates to be the same
  // equiv function of residue that two are close enough that they are equal
  // (instead of checking that UUID is the same)
  // Maybe this function exists already? Residue#is_equal, UUID false. May just work
  // ^^^ This is for each segment, they are in the structure object.
  // Step 1 is to write function that takes two segment graphs are equal within a threshold
  // Check numerical drift from rebuilding graph
  //
  // Weird edge cases
  // Need to put connections in database. This is a concept in segment_graph
  // Normally it's a parent/child connection, but start and end come back. To describe
  // that, there's a separate concept where two things are connected. There's an
  // add_connection method (interface?). There may be extra connections to recover.
  // When I'm writing it to the database, going to have to check if I've seen every
  // connection. Sorta like a sibling but not really. That connection needs
  // to be supported. Simple way of finding them is find connection where
  // neither are parents. Saved and readded. Absolutely critical for this to work.
  //
  // What to do with ref_frame and coordinate. Just for the roots,
  // if it's a root, translate it to that origin (whatever the coordinates are)
  // So when we load in (1), x, y, and z need to move to x_final, blah blah, move to saved origin
  // Have to use segment.move to do that is the difference between those x_initial to x_final (saved db) <--
  // For ref_frame, this is orientation of segment, defined by x, y, and z VECTORS.
  // Need to get from x_initial to x_final ... that are all vectors. Those vectors are stored as ref_frame
  // The rotational difference is a method called get_rotation_between_frames. I think this is segment.rotate or something
  // Rotation and then translation.
  //
  // Also want to move and rotate before saving
  // Want to add connections (I can just make shit up)
  // Add many different roots, things that aren't going to be added to the next one.
  // Just don't supply a connection. Think of edge cases.
  //
  // Ox DNA, this is 3D structure.

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
