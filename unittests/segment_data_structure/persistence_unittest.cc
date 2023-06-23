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

float compare_segment_coordinates(const Segment &s1, const Segment &s2) {
  const auto &structure1 = s1.get_structure();
  const auto &structure2 = s2.get_structure();
  std::vector<std::vector<std::string>> s1_coords;
  std::vector<std::vector<std::string>> s2_coords;

  // Was going to use this as a helper method to compare the differences
  // in coordinates for two segments. Delete this or work off of it.
  //
  // for (auto const &c : &s1.get_chains()) {
  //   auto chain_str = ::base::string::split(c.get_str(), ",");
  //   for (auto chain_data : chain_str) {
  //     auto chain_datum = ::base::string::split(chain_data, " ");
  //   }
  // }

  // for (auto const &c : &s2.get_chains()) {
  //   auto chain_str = ::base::string::split(c.get_str(), ",");
  //   for (auto chain_data : chain_str) {
  //     auto chain_datum = ::base::string::split(chain_data, " ");
  //   }
  // }
  return 0.0f;
}

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
      std::cout << "Is new graph the same as itself?\n";
      CHECK(db_graph == db_graph);
      std::cout << "Is new graph same when pulled from database again?\n";
      const segment_data_structure::SegmentGraphAllAtom second_db_graph = persistence.retrieve_segment_graph_from_database("test_seg", "test_dir/test.db");
      CHECK(db_graph == second_db_graph);
      std::cout << "Is in-memory graph same as db graph?\n";
      CHECK(sg == db_graph);
    }

    SUBCASE("Test coordinate drift") {
      const auto mem_seg1 = rm.get_segment(resource_management::SegmentInfo{"HELIX.IDEAL.6"});
      persistence.save_segment_to_database(
        *mem_seg1, "test_dir", "test", 999, "DatabaseSegment1"
      );
      const auto mem_seg2 = persistence.retrieve_segment_from_database("DatabaseSegment1", "test_dir/test.db");
      // Compare the coordinates
      // The plan here was to check the difference between each coordinate and then
      // to save mem_seg2 into the database, pull it out into a mem_seg3, and check
      // the difference there. Was going to do that repeatedly to see how much the
      // coordinates change each time they're saved to the database
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
}
