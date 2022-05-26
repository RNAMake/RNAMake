

#include "../common.hpp"
#include <doctest.h>
#include "math/xyz_matrix.hpp"
//#include "../base/file_io.h"
//#include "../base/settings.h"
//#include "../math/numerical.h"

/*
TEST_CASE("Test Matrix math ") {
  SUBCASE("Test Stringify Matrices") {
    auto m = math::Matrix3x3(5.0);
    auto s = m.get_str();
    auto m2 = math::Matrix3x3(s);

    CHECK(math::are_matrices_equal(m, m2));
  }

  SUBCASE("Single known test of unitarize compared to python") {
    auto path = base::unittest_resource_dir() + "/math/test_unitarize.dat";
    auto lines = base::get_lines_from_file(path);
    auto org_m = math::Matrix3x3(lines[0]);

    auto m = math::Matrix3x3(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 1.0, 1.0);

    auto unit = m.get_unitarize();

    CHECK(math::are_matrices_equal(org_m, unit));
  }

  SUBCASE("Test unitarize in batch with 1000 matrices") {
    auto path =
        base::unittest_resource_dir() + "/math/test_unitarize_multi.dat";
    auto lines = base::get_lines_from_file(path);

    int fail = 0;
    for (auto const& l : lines) {
      if (l.length() < 10) {
        break;
      }
      auto spl = base::split_str_by_delimiter(l, "|");
      auto org_m = math::Matrix3x3(spl[0]);
      auto final_m = math::Matrix3x3(spl[1]);

      auto unit = org_m.get_unitarize();

      if (!math::are_matrices_equal(final_m, unit)) {
        fail = 1;
        break;
      }

      org_m.unitarize();

      if (!math::are_matrices_equal(final_m, org_m)) {
        fail = 1;
        break;
      }
    }

    CHECK(fail == 0);
  }
}

*/