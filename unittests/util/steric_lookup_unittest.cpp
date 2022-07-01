#include "../common.hpp"

#include <util/point_generator.h>
#include <util/steric_lookup.hpp>

TEST_CASE("Test Steric Lookup for quick Sterics ") {
  auto p_generator = util::PointGenerator();
  SUBCASE("Test adding points to lookup") {
    auto sl = util::StericLookupNew();
    auto p = math::Vector3();
    sl.add_point(p);
    sl.to_pdb("grid.pdb");

    auto clash = 0, sl_clash = 0;
    double dist = 0;
    auto test_p = math::Vector3();

    // counter
    int count = 0;
    for (int i = 0; i < 1000; i++) {
      clash = 0;
      test_p = p_generator.rand_point(10);
      dist = p.distance(test_p);

      if (dist < 2.65) {
        clash = 1;
      }
      sl_clash = sl.clash(test_p);

      if (sl_clash == 0 && clash == 1) {
        count += 1;
      }
    }

    CHECK(count < 100);
  }
  SUBCASE("Test adding a set of points to lookup") {
    auto points = math::Vector3s();
    for (int i = 0; i < 1000; i++) {
      points.push_back(p_generator.rand_point(100));
    }

    auto sl = util::StericLookup();
    sl.add_points(points);

    auto dist = 0.0f;
    auto test_p = math::Vector3();
    auto clash = 0, miss_count = 0, sl_clash = 0;

    for (int i = 0; i < 10000; i++) {
      clash = 0;
      test_p = p_generator.rand_point(100);

      for (auto const &p : points) {
        dist = test_p.distance(p);
        if (dist < 2.65) {
          clash = 1;
          break;
        }
      }

      sl_clash = sl.clash(test_p);

      if (sl_clash != clash) {
        miss_count += 1;
      }
    }

    CHECK(miss_count < 200);
  }
  SUBCASE("test clash") {
    SUBCASE("test points directly on top of each other") {
      // add points to the lookup
      auto steric_lookup_test = util::StericLookupNew();
      auto point_1 = math::Vector3(1, 1, 1);
      steric_lookup_test.add_point(point_1);
      steric_lookup_test.to_pdb("test_clash.pdb");
      // plug them into the clash function
      auto point_2 = math::Vector3(1, 1, 1);
      bool do_points_clash = steric_lookup_test.clash(point_2);
      // check to see if clash spits out the right answer; true
      CHECK(do_points_clash == true);
    }
    SUBCASE("test points far away from each other") {
      // add points to the lookup
      auto steric_lookup_test = util::StericLookupNew();
      auto point_1 = math::Vector3(1, 1, 1);
      steric_lookup_test.add_point(point_1);
      steric_lookup_test.to_pdb("test_clash.pdb");
      // plug them into the clash function
      auto point_2 = math::Vector3(10, -3, 0);
      bool do_points_clash = steric_lookup_test.clash(point_2);
      // check to see if clash spits out the right answer; false
      CHECK(do_points_clash == false);
    }
    SUBCASE("test points that are really close to each other") {
      // TODO what is the acceptable radius?
      auto steric_lookup_test = util::StericLookupNew();
      auto point_1 = math::Vector3(1, 1, 1);
      steric_lookup_test.add_point(point_1);
      steric_lookup_test.to_pdb("test_clash.pdb");
      auto point_2 = math::Vector3(1.00000001, 1.00000001, 1.00000001);
      bool do_points_clash = steric_lookup_test.clash(point_2);
      CHECK(do_points_clash == true);
    }
    SUBCASE("test points that are really close to each other but far enough") {
      // TODO what is the acceptable radius?
      auto steric_lookup_test = util::StericLookupNew();
      auto point_1 = math::Vector3(1, 1, 1);
      steric_lookup_test.add_point(point_1);
      steric_lookup_test.to_pdb("test_clash.pdb");
      auto point_2 = math::Vector3(3, 3, 3);
      bool do_points_clash = steric_lookup_test.clash(point_2);
      CHECK(do_points_clash == false);
    }
  }
  SUBCASE("test clash between a set of points") {
    // simple formula for these tests:
    // initialize a list of points
    // add them into the lookup
    // plug them into the clash function
    SUBCASE("test clash between same points") {
      auto point_1 = math::Vector3(1, 1, 1);
      auto point_2 = math::Vector3(1, 1, 1);
      auto point_3 = math::Vector3(1, 1, 1);
      auto point_set_1 = math::Vector3s{point_1, point_2, point_3};

    }
    SUBCASE("test clash between completely different points") {

    }
    SUBCASE("test points that are really close to each other") {

    }
    SUBCASE("test points that are really close to each other but far enough") {

    }
  }
}
