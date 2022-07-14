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

    auto sl = util::StericLookupNew();
    sl.add_points(points);

    auto dist = 0.0f;
    auto test_p = math::Vector3();
    auto clash = 0, miss_count = 0, sl_clash = 0;

    for (int i = 0; i < 10000; i++) {
      clash = 0;
      test_p = p_generator.rand_point(100);

      for (auto const &p : points) {
        dist = test_p.distance(p);
        if (dist < 2.7) {
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
  SUBCASE("test clash points") {
    SUBCASE("test lone points") {
      SUBCASE("test points directly on top of each other") {
        auto steric_lookup_test = util::StericLookupNew();
        auto point_1 = math::Vector3(0, 0, 0);
        steric_lookup_test.add_point(point_1);
        auto point_2 = math::Vector3(0, 0, 0);
        CHECK(steric_lookup_test.clash(point_2) == true);
      }
      SUBCASE("test points far away from each other") {
        auto steric_lookup_test = util::StericLookupNew();
        auto point_1 = math::Vector3(0, 0, 0);
        steric_lookup_test.add_point(point_1);
        auto point_2 = math::Vector3(5, 0, 0);
        CHECK(steric_lookup_test.clash(point_2) == false);
      }
      SUBCASE("test points almost out of each other") {
        auto steric_lookup_test = util::StericLookupNew();
        auto point_1 = math::Vector3(0, 0, 0);
        steric_lookup_test.add_point(point_1);
        auto point_2 = math::Vector3(2.69, 0, 0);
        CHECK(steric_lookup_test.clash(point_2) == true);
      }
      SUBCASE("test points really close but far enough") {
        auto steric_lookup_test = util::StericLookupNew();
        auto point_1 = math::Vector3(0, 0, 0);
        steric_lookup_test.add_point(point_1);
        auto point_2 = math::Vector3(2.75, 0, 0);
        CHECK(steric_lookup_test.clash(point_2) == false);
      }
    }
    SUBCASE("test adding collection of points in xy plane") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(0, -1, 0);
      auto point_3 = math::Vector3(0, -2, 0);
      auto point_4 = math::Vector3(0, -3, 0);

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_point(point_1);
      steric_lookup_test.add_point(point_2);
      steric_lookup_test.add_point(point_3);
      steric_lookup_test.add_point(point_4);

      auto test_point_1 = math::Vector3(2.75, 0, 0);
      auto test_point_2 = math::Vector3(2.70, 0, 0);
      auto test_point_3 = math::Vector3(1.25, -1.25, 0);
      auto test_point_4 = math::Vector3(-1.75, 2, 0);
      auto test_point_5 = math::Vector3(-1.75, 2.1, 0);
      auto test_point_6 = math::Vector3(2.6, -1.5, 0);
      auto test_point_7 = math::Vector3(0, 0, 0);
      auto test_point_8 = math::Vector3(0, -1, 0);
      auto test_point_9 = math::Vector3(0, -2, 0);
      auto test_point_10 = math::Vector3(0, -3, 0);
      auto test_point_11 = math::Vector3(-2, 2, 0);

      CHECK(steric_lookup_test.clash(test_point_1) == false);
      CHECK(steric_lookup_test.clash(test_point_2) == true);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
      CHECK(steric_lookup_test.clash(test_point_4) == true);
      CHECK(steric_lookup_test.clash(test_point_5) == true);
      CHECK(steric_lookup_test.clash(test_point_6) == true);
      CHECK(steric_lookup_test.clash(test_point_7) == true);
      CHECK(steric_lookup_test.clash(test_point_8) == true);
      CHECK(steric_lookup_test.clash(test_point_9) == true);
      CHECK(steric_lookup_test.clash(test_point_10) == true);
      CHECK(steric_lookup_test.clash(test_point_11) == false);
    }
    SUBCASE("test adding collection of points in xz-plane") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(0, 0, 1);
      auto point_3 = math::Vector3(1, 0, 0);
      auto point_4 = math::Vector3(1, 0, 1);

      auto point_set_1 = math::Vector3s{point_1, point_2, point_3, point_4};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(0, 0, 1);
      auto test_point_3 = math::Vector3(1, 0, 0);
      auto test_point_4 = math::Vector3(1, 0, 1);
      auto test_point_5 = math::Vector3(3.75, 0, 0);
      auto test_point_6 = math::Vector3(0, 0, 3.75);
      auto test_point_7 = math::Vector3(3.7, 0, 0);
      auto test_point_8 = math::Vector3(0, 0, 3.7);
      auto test_point_9 = math::Vector3(1.25, 0, -1.25);
      auto test_point_10 = math::Vector3(0, 0, -2.75);
      auto test_point_11 = math::Vector3(0, 0, -2.69);
      auto test_point_12 = math::Vector3(0, 0, -2.75);

      CHECK(steric_lookup_test.clash(test_point_1) == true);
      CHECK(steric_lookup_test.clash(test_point_2) == true);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
      CHECK(steric_lookup_test.clash(test_point_4) == true);
      CHECK(steric_lookup_test.clash(test_point_5) == false);
      CHECK(steric_lookup_test.clash(test_point_6) == false);
      CHECK(steric_lookup_test.clash(test_point_7) == true);
      CHECK(steric_lookup_test.clash(test_point_8) == true);
      CHECK(steric_lookup_test.clash(test_point_9) == true);
      CHECK(steric_lookup_test.clash(test_point_10) == false);
      CHECK(steric_lookup_test.clash(test_point_11) == false);
      CHECK(steric_lookup_test.clash(test_point_12) == false);
    }
    SUBCASE("test adding collection of points in yz-plane") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(0, 0, 1);
      auto point_3 = math::Vector3(0, 1, 0);
      auto point_4 = math::Vector3(0, 1, 1);

      auto point_set_1 = math::Vector3s{point_1, point_2, point_3, point_4};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(0, 0, 1);
      auto test_point_3 = math::Vector3(0, 1, 0);
      auto test_point_4 = math::Vector3(0, 1, 1);
      auto test_point_5 = math::Vector3(0, 0, 3.75);

      CHECK(steric_lookup_test.clash(test_point_1) == true);
      CHECK(steric_lookup_test.clash(test_point_2) == true);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
      CHECK(steric_lookup_test.clash(test_point_4) == true);
      CHECK(steric_lookup_test.clash(test_point_5) == false);
    }
    SUBCASE("test adding collection of points in 3d space") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(1, 1, 1);
      auto point_3 = math::Vector3(2, 2, 2);
      auto point_4 = math::Vector3(3, 3, 3);

      auto point_set_1 = math::Vector3s{point_1, point_2, point_3, point_4};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(1, 1, 1);
      auto test_point_3 = math::Vector3(2, 2, 2);
      auto test_point_4 = math::Vector3(3, 3, 3);
      auto test_point_5 = math::Vector3(-1.75, -1.75, -1.75);
      auto test_point_6 = math::Vector3(-1.5, -1.5, -1.5);
      auto test_point_7 = math::Vector3(-1.6, -1.6, -1.6);
      auto test_point_8 = math::Vector3(4.25, 4.25, 4.25);
      auto test_point_9 = math::Vector3(0, 0, 3.25);
      auto test_point_10 = math::Vector3(0, 0.5, 3.25);
      auto test_point_11 = math::Vector3(0, 0, 3.5);

      CHECK(steric_lookup_test.clash(test_point_1) == true);
      CHECK(steric_lookup_test.clash(test_point_2) == true);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
      CHECK(steric_lookup_test.clash(test_point_4) == true);
      CHECK(steric_lookup_test.clash(test_point_5) == false);
      CHECK(steric_lookup_test.clash(test_point_6) == true);
      CHECK(steric_lookup_test.clash(test_point_7) == false);
      CHECK(steric_lookup_test.clash(test_point_8) == true);
      CHECK(steric_lookup_test.clash(test_point_9) == true);
      CHECK(steric_lookup_test.clash(test_point_10) == true);
      CHECK(steric_lookup_test.clash(test_point_11) == false);
    }
  }
  SUBCASE("test clash between a set of points") {
    SUBCASE("test clash between same points") {
      auto point_1 = math::Vector3(1, 1, 1);
      auto point_2 = math::Vector3(1, 1, 1);
      auto point_3 = math::Vector3(1, 1, 1);
      auto point_set_1 = math::Vector3s{point_1, point_2, point_3};
      auto point_set_2 = math::Vector3s{point_1, point_2};
      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);
      CHECK(steric_lookup_test.clash(point_set_2) == true);
    }
    SUBCASE("test clash between completely different points") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(0, -1, 0);
      auto point_3 = math::Vector3(0, -2, 0);
      auto point_set_1 = math::Vector3s{point_1, point_2, point_3};
      auto point_4 = math::Vector3(3, 1, 0);
      auto point_5 = math::Vector3(-4, 2, -4);
      auto point_set_2 = math::Vector3s{point_4, point_5};
      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);
      CHECK(steric_lookup_test.clash(point_set_2) == false);
    }
    SUBCASE("test points that are really close to each other") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(1, 0, 0);
      auto point_3 = math::Vector3(2, 0, 0);
      auto point_4 = math::Vector3(3, 0, 0);
      auto point_set_1 = math::Vector3s{point_1, point_2, point_3, point_4};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto point_5 = math::Vector3(0, -3, 2);
      auto point_6 = math::Vector3(1, -1, 0);
      auto point_7 = math::Vector3(0, 1, 1);
      auto point_8 = math::Vector3(1.3, 0, -1.2);

      auto point_9 = math::Vector3(5, 0, 0);
      auto point_10 = math::Vector3(2, -2.7, 0);
      auto point_11 = math::Vector3(3.1, 0, -2.4);

      auto point_set_2 = math::Vector3s{point_5, point_6, point_7, point_8};
      auto point_set_3 = math::Vector3s{point_9, point_10};

      CHECK(steric_lookup_test.clash(point_set_2) == true); // should be true
      CHECK(steric_lookup_test.clash(point_11) == true);
    }
    SUBCASE("test points that are really close to boundary but far enough") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_2 = math::Vector3(1, 0, 0);
      auto point_3 = math::Vector3(2, 0, 0);
      auto point_4 = math::Vector3(3, 0, 0);
      auto point_set_1 = math::Vector3s{point_1, point_2, point_3, point_4};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto point_5 = math::Vector3(6, 0, 0);
      auto point_6 = math::Vector3(5.75, 0, 0);
      auto point_7 = math::Vector3(-2.75, 0, 0);

      auto point_set_2 = math::Vector3s{point_5, point_6, point_7};

      auto point_8 = math::Vector3(3, -3, 0);
      auto point_9 = math::Vector3(3, -2.75, 0);
      auto point_10 = math::Vector3(-2.75, 0, 0);

      auto point_set_3 = math::Vector3s{point_8, point_9, point_10};

      CHECK(steric_lookup_test.clash(point_set_2) == false);
      CHECK(steric_lookup_test.clash(point_set_3) == false);
    }
    SUBCASE("test mixture of points in and outside radius") {
      auto point_1 = math::Vector3(0, 0, 0);
      auto point_set_1 = math::Vector3s{point_1};

      auto steric_lookup_test = util::StericLookupNew();
      steric_lookup_test.add_points(point_set_1);

      auto point_5 = math::Vector3(5.75, 0, 0);       // outside
      auto point_6 = math::Vector3(1.75, 1.75, 1.75); // outside
      auto point_7 = math::Vector3(-3, -3, -1);       // outside

      auto point_8 = math::Vector3(0.21, -0.1, -0.001); // inside
      auto point_9 = math::Vector3(-0.52, 0.4, 0.1);    // inside
      auto point_10 = math::Vector3(2, 1, 0);           // inside
      auto point_11 = math::Vector3(1.55, 1.55, 1.55);  // inside

      auto point_12 = math::Vector3(1, 1, 1);
      auto point_13 = math::Vector3(0.75, 0, -1);

      auto point_set_2 = math::Vector3s{point_5, point_6, point_7};
      auto point_set_3 = math::Vector3s{point_8, point_9, point_10, point_11};

      auto point_set_4 = math::Vector3s{point_5,  point_8,  point_7,
                                        point_11, point_12, point_13};

      CHECK(steric_lookup_test.clash(point_set_2) == false);
      CHECK(steric_lookup_test.clash(point_set_3) == true);
      CHECK(steric_lookup_test.clash(point_set_4) == true);
    }
  }
  SUBCASE("test steric lookups around zero") {
    SUBCASE("test lookups with point 0, 0, 0") {
      auto point_1 = math::Vector3(0, 0, 0);

      auto test_steric_zero = util::StericLookupNew();
      test_steric_zero.add_point(point_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(0.15, 0.15, 0);
      auto test_point_3 = math::Vector3(0.4, 0, 0);
      auto test_point_4 = math::Vector3(0.6, -0.3, -0.1);
      auto test_point_5 = math::Vector3(-0.75, 0, -1);

      CHECK(test_steric_zero.clash(test_point_1) == true);
      CHECK(test_steric_zero.clash(test_point_2) == true);
      CHECK(test_steric_zero.clash(test_point_3) == true);
      CHECK(test_steric_zero.clash(test_point_4) == true);
      CHECK(test_steric_zero.clash(test_point_5) == true);
    }
    SUBCASE("test lookups with point boundary near origin") {
      auto point_1 = math::Vector3(3.5, 0, 0);

      auto test_steric_zero = util::StericLookupNew();
      test_steric_zero.add_point(point_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(0.75, 0, 0);
      auto test_point_3 = math::Vector3(1.1, 0, 0);
      auto test_point_4 = math::Vector3(0.21, -0.3, 1.1);
      auto test_point_5 = math::Vector3(-0.5, 0.01, 0.15);

      CHECK(test_steric_zero.clash(test_point_1) == false);
      CHECK(test_steric_zero.clash(test_point_2) == false);
      CHECK(test_steric_zero.clash(test_point_3) == true);
      CHECK(test_steric_zero.clash(test_point_4) == false);
      CHECK(test_steric_zero.clash(test_point_5) == false);
    }
  }
  SUBCASE("test lookups near the boundaries") {
    SUBCASE("test lower boundary") {
      SUBCASE("test lookup boundary grazing the bounding box") {
        auto point_1 = math::Vector3(-197.25, -197.25, -197.25);

        auto steric_lookup_test = util::StericLookupNew();
        steric_lookup_test.add_point(point_1);

        auto test_point_1 = math::Vector3(-198, -198, -198);
        auto test_point_2 = math::Vector3(-199, -199, -199);
        auto test_point_3 = math::Vector3(-200, -200, -200);

        CHECK(steric_lookup_test.clash(test_point_1) == true);
        CHECK(steric_lookup_test.clash(test_point_2) == false);
        CHECK(steric_lookup_test.clash(test_point_3) == false);
      }
      SUBCASE("test lookups completely outside the bounding box") {
        auto point_1 = math::Vector3(-200.25, -200.25, -200.25);
        auto steric_lookup_test = util::StericLookupNew();
        CHECK_THROWS_AS(steric_lookup_test.add_point(point_1),
                        base::MathException);
      }
    }
    SUBCASE("test upper boundary cases") {
      SUBCASE("test lookup boundary grazing the bounding box") {
        auto point_1 = math::Vector3(99.75, 99.75, 99.75);

        auto steric_lookup_test = util::StericLookupNew();
        steric_lookup_test.add_point(point_1);

        auto test_point_1 = math::Vector3(99, 99, 99);
        auto test_point_2 = math::Vector3(100, 100, 100);
        auto test_point_3 = math::Vector3(101, 101, 101);

        CHECK(steric_lookup_test.clash(test_point_1) == true);
        CHECK(steric_lookup_test.clash(test_point_2) == true);
        CHECK(steric_lookup_test.clash(test_point_3) == true);
      }
      SUBCASE("test lookup completely outside the bounding box") {
        auto point_1 = math::Vector3(100.25, 100.25, 100.25);
        auto steric_lookup_test = util::StericLookupNew();
        CHECK_THROWS_AS(steric_lookup_test.add_point(point_1),
                        base::MathException);
      }
    }
  }

  SUBCASE("test custom steric lookups") {
    SUBCASE("test lookup with radius of 3") {
      auto point_1 = math::Vector3(0, 0, 0);

      // grid_size 0.5, cutoff 3, radius 12
      auto steric_lookup_test = util::StericLookupNew(0.5, 3, 12);
      steric_lookup_test.add_point(point_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(3.0, 0, 0);
      auto test_point_3 = math::Vector3(2.75, 0, 0);
      auto test_point_4 = math::Vector3(3.01, 0, 0);
      auto test_point_5 = math::Vector3(2.99, 0, 0);
      auto test_point_6 = math::Vector3(2.9995, 0, 0);

      CHECK(steric_lookup_test.clash(test_point_1) == true);
      // TODO check boundaries greater/less than
      // CHECK(steric_lookup_test.clash(test_point_2) == false);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
      CHECK(steric_lookup_test.clash(test_point_4) == false);
      CHECK(steric_lookup_test.clash(test_point_5) == true);
      CHECK(steric_lookup_test.clash(test_point_6) == true);
    }
    SUBCASE("test lookup with radius of 1") {
      auto point_1 = math::Vector3(0, 0, 0);

      // grid_size 0.1, cutoff 1, radius 12
      auto steric_lookup_test = util::StericLookupNew(0.1, 1, 12);
      steric_lookup_test.add_point(point_1);

      auto test_point_1 = math::Vector3(0, 0, 0);
      auto test_point_2 = math::Vector3(0.999, 0, 0);
      auto test_point_3 = math::Vector3(0.577, 0.577, 0.577);

      CHECK(steric_lookup_test.clash(test_point_1) == true);
      CHECK(steric_lookup_test.clash(test_point_2) == true);
      CHECK(steric_lookup_test.clash(test_point_3) == true);
    }
    SUBCASE("test negative arguments in the lookup") {
      CHECK_THROWS_AS(util::StericLookupNew(-0.5, 2, 12), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(0.3, -1, 12), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(0.25, 4, -12), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(-0.2, -10, -12),
                      base::MathException);
    }
    SUBCASE("test zero arguments in lookup") {
      CHECK_THROWS_AS(util::StericLookupNew(0, 2, 12), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(0.3, 0, 12), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(0.25, 4, 0), base::MathException);
      CHECK_THROWS_AS(util::StericLookupNew(0, 0, 0), base::MathException);
    }
  }
  SUBCASE("test counting of clashes") {
    auto point_1 = math::Vector3(0, 0, 0);

    auto steric_lookup_test = util::StericLookupNew();
    steric_lookup_test.add_point(point_1);

    auto test_point_1 = math::Vector3(0, 0, 0);            // true
    auto test_point_2 = math::Vector3(1, 1, 1);            // true
    auto test_point_3 = math::Vector3(2, 2, 2);            // false
    auto test_point_4 = math::Vector3(5, 5, 5);            // false
    auto test_point_5 = math::Vector3(-4, 1, 3);           // false
    auto test_point_6 = math::Vector3(0.21, 1, 0);         // true
    auto test_point_7 = math::Vector3(-1.11, 0.22, 0.45);  // true
    auto test_point_8 = math::Vector3(-0.34, -1, 0);       // true
    auto test_point_9 = math::Vector3(0, -1, 0);           // true
    auto test_point_10 = math::Vector3(-1.75, 1.75, 1.75); // false

    auto point_set_1 = math::Vector3s{
        test_point_1, test_point_2, test_point_3, test_point_4, test_point_5,
        test_point_6, test_point_7, test_point_8, test_point_9, test_point_10};

    CHECK(steric_lookup_test.clash(test_point_1) == true);
    CHECK(steric_lookup_test.clash(test_point_2) == true);
    CHECK(steric_lookup_test.clash(test_point_3) == false);
    CHECK(steric_lookup_test.clash(test_point_4) == false);
    CHECK(steric_lookup_test.clash(test_point_5) == false);
    CHECK(steric_lookup_test.clash(test_point_6) == true);
    CHECK(steric_lookup_test.clash(test_point_7) == true);
    CHECK(steric_lookup_test.clash(test_point_8) == true);
    CHECK(steric_lookup_test.clash(test_point_9) == true);
    CHECK(steric_lookup_test.clash(test_point_10) == false);

    CHECK(steric_lookup_test.total_clash(point_set_1) == doctest::Approx(6));
  }

  SUBCASE("") {

  }

}

// TODO think of unittests