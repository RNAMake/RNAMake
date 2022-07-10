//
// Created by Joseph Yesselman on 1/3/18.
//

#include "../common.hpp"
#include <iostream>

#include <base/paths.hpp>
#include <resource_management/segment_sqlite_library.h>
#include <resource_management/sqlite_library.h>

using namespace resource_management;

TEST_CASE("Test basic sqlite library") {
  SUBCASE("test base") {
    auto db_path =
        base::path::unittest_resource_path() + "/resource_management/test.db";
    auto sqlib = SqliteLibrary(db_path, "data_table");
    SUBCASE("database file must exist") {
      REQUIRE_THROWS_AS(SqliteLibrary("test_fake.db", "data_table"),
                        ResourceManagementException);
    }
  }

  SUBCASE("test segment databases") {
    auto db_path =
        base::path::resources_path() + "/motif_libraries_new/ideal_helices.db";
    auto seg_lib = SegmentSqliteLibrary(db_path, "data_table");
    auto seg = seg_lib.get_segment(SegmentInfo{"HELIX.IDEAL.2"});
    REQUIRE(seg.get_name() == "HELIX.IDEAL.2");

    auto seg3 = seg_lib.get_segment(SegmentInfo{"", "", "GGGG_LLLL_CCCC_RRRR"});

    // REQUIRE(seg->is_equal(*seg3, false));
  }

  /*SECTION("test retrivial is the same from just using seg factory") {
    auto seg_factory = all_atom::SegmentFactory(rts);
    auto seg =
        seg_lib.get_segment(StringStringMap{{"name", "HELIX.IDEAL.2"}});
    auto path = base::unittest_resources_path() +
                "/all_atom/HELIX.IDEAL.2/HELIX.IDEAL.2.pdb";
    auto seg2 =
        seg_factory.segment_from_pdb(path, util::SegmentType::HELIX, true);
    seg_factory.align_segment_to_ref_frame(*seg2);

    REQUIRE(seg->is_equal(*seg2, false));
  }  */

  SUBCASE("test errors for invalid queries") {
    auto db_path =
        base::path::resources_path() + "/motif_libraries_new/ideal_helices.db";
    auto seg_lib = SegmentSqliteLibrary(db_path, "data_table");
    REQUIRE_THROWS_AS(seg_lib.get_segment(SegmentInfo{"FAKE"}),
                      ResourceManagementException);

    /*REQUIRE_THROWS_AS(seg_lib.get_segment(StringStringMap{{"end_id",
"FAKE"}}), resource_management::SqliteLibraryException);

    REQUIRE_THROWS_AS(seg_lib.get_segment(StringStringMap{{"end_name",
"FAKE"}}), resource_management::SqliteLibraryException);

    REQUIRE_THROWS_AS(seg_lib.get_segment(StringStringMap{{"name",
"HELIX.IDEAL.2"}, {"end_id", "FAKE"}}),
                      resource_management::SqliteLibraryException); */
  }

  SUBCASE("test contains segments") {
    auto db_path =
        base::path::resources_path() + "/motif_libraries_new/ideal_helices.db";
    auto seg_lib = SegmentSqliteLibrary(db_path, "data_table");
    REQUIRE(seg_lib.contains_segment(SegmentInfo{"HELIX.IDEAL.2"}) == true);
    REQUIRE(seg_lib.contains_segment(SegmentInfo{"FAKE"}) == false);
  }


}