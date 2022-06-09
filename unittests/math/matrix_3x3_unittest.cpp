#include "../common.hpp"
#include <sstream>
#include <math/matrix_3x3.hpp>
//#include <doctest.h>

TEST_CASE("Test Matrix math ") {
  /*
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


    CHECK(fail == 0);
  }
}

*/
  SUBCASE("test constructors") {
    SUBCASE("test empty") {
      math::Matrix3x3 matrix_1 = {};
      CHECK(matrix_1.get_xx() == doctest::Approx(0));
      CHECK(matrix_1.get_xy() == doctest::Approx(0));
      CHECK(matrix_1.get_xz() == doctest::Approx(0));
      CHECK(matrix_1.get_yx() == doctest::Approx(0));
      CHECK(matrix_1.get_yy() == doctest::Approx(0));
      CHECK(matrix_1.get_yz() == doctest::Approx(0));
      CHECK(matrix_1.get_zx() == doctest::Approx(0));
      CHECK(matrix_1.get_zy() == doctest::Approx(0));
      CHECK(matrix_1.get_zz() == doctest::Approx(0));
    }
    SUBCASE("test zero") {
      math::Matrix3x3 matrix_1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      CHECK(matrix_1.get_xx() == doctest::Approx(0));
      CHECK(matrix_1.get_xy() == doctest::Approx(0));
      CHECK(matrix_1.get_xz() == doctest::Approx(0));
      CHECK(matrix_1.get_yx() == doctest::Approx(0));
      CHECK(matrix_1.get_yy() == doctest::Approx(0));
      CHECK(matrix_1.get_yz() == doctest::Approx(0));
      CHECK(matrix_1.get_zx() == doctest::Approx(0));
      CHECK(matrix_1.get_zy() == doctest::Approx(0));
      CHECK(matrix_1.get_zz() == doctest::Approx(0));
    }
    SUBCASE("test simple matrix") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      CHECK(matrix_1.get_xx() == doctest::Approx(1));
      CHECK(matrix_1.get_xy() == doctest::Approx(2));
      CHECK(matrix_1.get_xz() == doctest::Approx(3));
      CHECK(matrix_1.get_yx() == doctest::Approx(4));
      CHECK(matrix_1.get_yy() == doctest::Approx(5));
      CHECK(matrix_1.get_yz() == doctest::Approx(6));
      CHECK(matrix_1.get_zx() == doctest::Approx(7));
      CHECK(matrix_1.get_zy() == doctest::Approx(8));
      CHECK(matrix_1.get_zz() == doctest::Approx(9));
    }
    SUBCASE("test complex matrix") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      CHECK(matrix_1.get_xx() == doctest::Approx(54));
      CHECK(matrix_1.get_xy() == doctest::Approx(23));
      CHECK(matrix_1.get_xz() == doctest::Approx(-33));
      CHECK(matrix_1.get_yx() == doctest::Approx(42));
      CHECK(matrix_1.get_yy() == doctest::Approx(-5));
      CHECK(matrix_1.get_yz() == doctest::Approx(16));
      CHECK(matrix_1.get_zx() == doctest::Approx(-17));
      CHECK(matrix_1.get_zy() == doctest::Approx(38));
      CHECK(matrix_1.get_zz() == doctest::Approx(19));
    }
    SUBCASE("test uniform value constructor") {
      math::Matrix3x3 matrix_1(1);
      CHECK(matrix_1.get_xx() == doctest::Approx(1));
      CHECK(matrix_1.get_xy() == doctest::Approx(1));
      CHECK(matrix_1.get_xz() == doctest::Approx(1));
      CHECK(matrix_1.get_yx() == doctest::Approx(1));
      CHECK(matrix_1.get_yy() == doctest::Approx(1));
      CHECK(matrix_1.get_yz() == doctest::Approx(1));
      CHECK(matrix_1.get_zx() == doctest::Approx(1));
      CHECK(matrix_1.get_zy() == doctest::Approx(1));
      CHECK(matrix_1.get_zz() == doctest::Approx(1));
    }
    SUBCASE("test complex uniform value constructior") {
      math::Matrix3x3 matrix_1(-0.22);
      CHECK(matrix_1.get_xx() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_xy() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_xz() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_yx() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_yy() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_yz() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_zx() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_zy() == doctest::Approx(-0.22));
      CHECK(matrix_1.get_zz() == doctest::Approx(-0.22));
    }
  }
  SUBCASE("test value/matrix addition") {
    SUBCASE("test zero value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 0;
      math::Matrix3x3 matrix_2 = t + matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(1));
      CHECK(matrix_2.get_xy() == doctest::Approx(2));
      CHECK(matrix_2.get_xz() == doctest::Approx(3));
      CHECK(matrix_2.get_yx() == doctest::Approx(4));
      CHECK(matrix_2.get_yy() == doctest::Approx(5));
      CHECK(matrix_2.get_yz() == doctest::Approx(6));
      CHECK(matrix_2.get_zx() == doctest::Approx(7));
      CHECK(matrix_2.get_zy() == doctest::Approx(8));
      CHECK(matrix_2.get_zz() == doctest::Approx(9));
    }
    SUBCASE("test zero matrix") {
      math::Matrix3x3 matrix_1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      Real t = 2;
      math::Matrix3x3 matrix_2 = t + matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(2));
      CHECK(matrix_2.get_xy() == doctest::Approx(2));
      CHECK(matrix_2.get_xz() == doctest::Approx(2));
      CHECK(matrix_2.get_yx() == doctest::Approx(2));
      CHECK(matrix_2.get_yy() == doctest::Approx(2));
      CHECK(matrix_2.get_yz() == doctest::Approx(2));
      CHECK(matrix_2.get_zx() == doctest::Approx(2));
      CHECK(matrix_2.get_zy() == doctest::Approx(2));
      CHECK(matrix_2.get_zz() == doctest::Approx(2));
    }
    SUBCASE("test simple matrix") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 1;
      math::Matrix3x3 matrix_2 = t + matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(2));
      CHECK(matrix_2.get_xy() == doctest::Approx(3));
      CHECK(matrix_2.get_xz() == doctest::Approx(4));
      CHECK(matrix_2.get_yx() == doctest::Approx(5));
      CHECK(matrix_2.get_yy() == doctest::Approx(6));
      CHECK(matrix_2.get_yz() == doctest::Approx(7));
      CHECK(matrix_2.get_zx() == doctest::Approx(8));
      CHECK(matrix_2.get_zy() == doctest::Approx(9));
      CHECK(matrix_2.get_zz() == doctest::Approx(10));
    }
    SUBCASE("test complex matrix") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = 10;
      math::Matrix3x3 matrix_2 = t + matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(64));
      CHECK(matrix_2.get_xy() == doctest::Approx(33));
      CHECK(matrix_2.get_xz() == doctest::Approx(-23));
      CHECK(matrix_2.get_yx() == doctest::Approx(52));
      CHECK(matrix_2.get_yy() == doctest::Approx(5));
      CHECK(matrix_2.get_yz() == doctest::Approx(26));
      CHECK(matrix_2.get_zx() == doctest::Approx(-7));
      CHECK(matrix_2.get_zy() == doctest::Approx(48));
      CHECK(matrix_2.get_zz() == doctest::Approx(29));
    }
  }
  SUBCASE("test matrix/value subtraction") {
    SUBCASE("test zero value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 0;
      math::Matrix3x3 matrix_2 = t - matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(-1));
      CHECK(matrix_2.get_xy() == doctest::Approx(-2));
      CHECK(matrix_2.get_xz() == doctest::Approx(-3));
      CHECK(matrix_2.get_yx() == doctest::Approx(-4));
      CHECK(matrix_2.get_yy() == doctest::Approx(-5));
      CHECK(matrix_2.get_yz() == doctest::Approx(-6));
      CHECK(matrix_2.get_zx() == doctest::Approx(-7));
      CHECK(matrix_2.get_zy() == doctest::Approx(-8));
      CHECK(matrix_2.get_zz() == doctest::Approx(-9));
    }
    SUBCASE("test zero matrix") {
      math::Matrix3x3 matrix_1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      Real t = 2;
      math::Matrix3x3 matrix_2 = matrix_1 - 2;
      CHECK(matrix_2.get_xx() == doctest::Approx(-2));
      CHECK(matrix_2.get_xy() == doctest::Approx(-2));
      CHECK(matrix_2.get_xz() == doctest::Approx(-2));
      CHECK(matrix_2.get_yx() == doctest::Approx(-2));
      CHECK(matrix_2.get_yy() == doctest::Approx(-2));
      CHECK(matrix_2.get_yz() == doctest::Approx(-2));
      CHECK(matrix_2.get_zx() == doctest::Approx(-2));
      CHECK(matrix_2.get_zy() == doctest::Approx(-2));
      CHECK(matrix_2.get_zz() == doctest::Approx(-2));
    }
    SUBCASE("test simple matrix") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 1;
      math::Matrix3x3 matrix_2 = matrix_1 - 1;
      CHECK(matrix_2.get_xx() == doctest::Approx(0));
      CHECK(matrix_2.get_xy() == doctest::Approx(1));
      CHECK(matrix_2.get_xz() == doctest::Approx(2));
      CHECK(matrix_2.get_yx() == doctest::Approx(3));
      CHECK(matrix_2.get_yy() == doctest::Approx(4));
      CHECK(matrix_2.get_yz() == doctest::Approx(5));
      CHECK(matrix_2.get_zx() == doctest::Approx(6));
      CHECK(matrix_2.get_zy() == doctest::Approx(7));
      CHECK(matrix_2.get_zz() == doctest::Approx(8));
    }
    SUBCASE("test complex matrix") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = 10;
      math::Matrix3x3 matrix_2 = matrix_1 - t;
      CHECK(matrix_2.get_xx() == doctest::Approx(44));
      CHECK(matrix_2.get_xy() == doctest::Approx(13));
      CHECK(matrix_2.get_xz() == doctest::Approx(-43));
      CHECK(matrix_2.get_yx() == doctest::Approx(32));
      CHECK(matrix_2.get_yy() == doctest::Approx(-15));
      CHECK(matrix_2.get_yz() == doctest::Approx(6));
      CHECK(matrix_2.get_zx() == doctest::Approx(-27));
      CHECK(matrix_2.get_zy() == doctest::Approx(28));
      CHECK(matrix_2.get_zz() == doctest::Approx(9));
    }
  }
  SUBCASE("test value/matrix multiplication") {
    SUBCASE("test zero value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 0;
      math::Matrix3x3 matrix_2 = t * matrix_1;
      CHECK(matrix_2.get_xx() == doctest::Approx(0));
      CHECK(matrix_2.get_xy() == doctest::Approx(0));
      CHECK(matrix_2.get_xz() == doctest::Approx(0));
      CHECK(matrix_2.get_yx() == doctest::Approx(0));
      CHECK(matrix_2.get_yy() == doctest::Approx(0));
      CHECK(matrix_2.get_yz() == doctest::Approx(0));
      CHECK(matrix_2.get_zx() == doctest::Approx(0));
      CHECK(matrix_2.get_zy() == doctest::Approx(0));
      CHECK(matrix_2.get_zz() == doctest::Approx(0));
    }
    SUBCASE("test zero matrix") {
      math::Matrix3x3 matrix_1 = {0, 0, 0, 0, 0, 0, 0, 0, 0};
      Real t = 2;
      math::Matrix3x3 matrix_2 = matrix_1 * 1;
      CHECK(matrix_2.get_xx() == doctest::Approx(0));
      CHECK(matrix_2.get_xy() == doctest::Approx(0));
      CHECK(matrix_2.get_xz() == doctest::Approx(0));
      CHECK(matrix_2.get_yx() == doctest::Approx(0));
      CHECK(matrix_2.get_yy() == doctest::Approx(0));
      CHECK(matrix_2.get_yz() == doctest::Approx(0));
      CHECK(matrix_2.get_zx() == doctest::Approx(0));
      CHECK(matrix_2.get_zy() == doctest::Approx(0));
      CHECK(matrix_2.get_zz() == doctest::Approx(0));
    }
    SUBCASE("test simple matrix/simple value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 1;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(1));
      CHECK(matrix_2.get_xy() == doctest::Approx(2));
      CHECK(matrix_2.get_xz() == doctest::Approx(3));
      CHECK(matrix_2.get_yx() == doctest::Approx(4));
      CHECK(matrix_2.get_yy() == doctest::Approx(5));
      CHECK(matrix_2.get_yz() == doctest::Approx(6));
      CHECK(matrix_2.get_zx() == doctest::Approx(7));
      CHECK(matrix_2.get_zy() == doctest::Approx(8));
      CHECK(matrix_2.get_zz() == doctest::Approx(9));
    }
    SUBCASE("test simple matrix/complex value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 10;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(10));
      CHECK(matrix_2.get_xy() == doctest::Approx(20));
      CHECK(matrix_2.get_xz() == doctest::Approx(30));
      CHECK(matrix_2.get_yx() == doctest::Approx(40));
      CHECK(matrix_2.get_yy() == doctest::Approx(50));
      CHECK(matrix_2.get_yz() == doctest::Approx(60));
      CHECK(matrix_2.get_zx() == doctest::Approx(70));
      CHECK(matrix_2.get_zy() == doctest::Approx(80));
      CHECK(matrix_2.get_zz() == doctest::Approx(90));
    }
    SUBCASE("test simple matrix/negative value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = -10;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(-10));
      CHECK(matrix_2.get_xy() == doctest::Approx(-20));
      CHECK(matrix_2.get_xz() == doctest::Approx(-30));
      CHECK(matrix_2.get_yx() == doctest::Approx(-40));
      CHECK(matrix_2.get_yy() == doctest::Approx(-50));
      CHECK(matrix_2.get_yz() == doctest::Approx(-60));
      CHECK(matrix_2.get_zx() == doctest::Approx(-70));
      CHECK(matrix_2.get_zy() == doctest::Approx(-80));
      CHECK(matrix_2.get_zz() == doctest::Approx(-90));
    }
    SUBCASE("test simple matrix/decimal value") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      Real t = 0.25;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(0.25));
      CHECK(matrix_2.get_xy() == doctest::Approx(0.5));
      CHECK(matrix_2.get_xz() == doctest::Approx(0.75));
      CHECK(matrix_2.get_yx() == doctest::Approx(1));
      CHECK(matrix_2.get_yy() == doctest::Approx(1.25));
      CHECK(matrix_2.get_yz() == doctest::Approx(1.5));
      CHECK(matrix_2.get_zx() == doctest::Approx(1.75));
      CHECK(matrix_2.get_zy() == doctest::Approx(2));
      CHECK(matrix_2.get_zz() == doctest::Approx(2.25));
    }
    SUBCASE("test complex matrix/simple value") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = 1;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(54));
      CHECK(matrix_2.get_xy() == doctest::Approx(23));
      CHECK(matrix_2.get_xz() == doctest::Approx(-33));
      CHECK(matrix_2.get_yx() == doctest::Approx(42));
      CHECK(matrix_2.get_yy() == doctest::Approx(-5));
      CHECK(matrix_2.get_yz() == doctest::Approx(16));
      CHECK(matrix_2.get_zx() == doctest::Approx(-17));
      CHECK(matrix_2.get_zy() == doctest::Approx(38));
      CHECK(matrix_2.get_zz() == doctest::Approx(19));
    }
    SUBCASE("test complex matrix/complex value") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = 17;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(918));
      CHECK(matrix_2.get_xy() == doctest::Approx(391));
      CHECK(matrix_2.get_xz() == doctest::Approx(-561));
      CHECK(matrix_2.get_yx() == doctest::Approx(714));
      CHECK(matrix_2.get_yy() == doctest::Approx(-85));
      CHECK(matrix_2.get_yz() == doctest::Approx(272));
      CHECK(matrix_2.get_zx() == doctest::Approx(-289));
      CHECK(matrix_2.get_zy() == doctest::Approx(646));
      CHECK(matrix_2.get_zz() == doctest::Approx(323));
    }
    SUBCASE("test complex matrix/negative value") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = -1;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(-54));
      CHECK(matrix_2.get_xy() == doctest::Approx(-23));
      CHECK(matrix_2.get_xz() == doctest::Approx(33));
      CHECK(matrix_2.get_yx() == doctest::Approx(-42));
      CHECK(matrix_2.get_yy() == doctest::Approx(5));
      CHECK(matrix_2.get_yz() == doctest::Approx(-16));
      CHECK(matrix_2.get_zx() == doctest::Approx(17));
      CHECK(matrix_2.get_zy() == doctest::Approx(-38));
      CHECK(matrix_2.get_zz() == doctest::Approx(-19));
    }
    SUBCASE("test complex matrix/decimal value") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      Real t = 0.25;
      math::Matrix3x3 matrix_2 = matrix_1 * t;
      CHECK(matrix_2.get_xx() == doctest::Approx(13.5));
      CHECK(matrix_2.get_xy() == doctest::Approx(5.75));
      CHECK(matrix_2.get_xz() == doctest::Approx(-8.25));
      CHECK(matrix_2.get_yx() == doctest::Approx(10.5));
      CHECK(matrix_2.get_yy() == doctest::Approx(-1.25));
      CHECK(matrix_2.get_yz() == doctest::Approx(4));
      CHECK(matrix_2.get_zx() == doctest::Approx(-4.25));
      CHECK(matrix_2.get_zy() == doctest::Approx(9.5));
      CHECK(matrix_2.get_zz() == doctest::Approx(4.75));
    }
  }
  SUBCASE("test division") {
    SUBCASE("test divide by zero") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      CHECK_THROWS_AS(matrix_1 / 0, base::MathException);
    }
    SUBCASE("test simple division") {
      math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 0};
      math::Matrix3x3 matrix_2 = matrix_1 / 10;
      CHECK(matrix_2.get_xx() == doctest::Approx(0.1));
      CHECK(matrix_2.get_xy() == doctest::Approx(0.2));
      CHECK(matrix_2.get_xz() == doctest::Approx(0.3));
      CHECK(matrix_2.get_yx() == doctest::Approx(0.4));
      CHECK(matrix_2.get_yy() == doctest::Approx(0.5));
      CHECK(matrix_2.get_yz() == doctest::Approx(0.6));
      CHECK(matrix_2.get_zx() == doctest::Approx(0.7));
      CHECK(matrix_2.get_zy() == doctest::Approx(0.8));
      CHECK(matrix_2.get_zz() == doctest::Approx(0));
    }
    SUBCASE("test complex division") {
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      math::Matrix3x3 matrix_2 = matrix_1 / 4;
      CHECK(matrix_2.get_xx() == doctest::Approx(13.5));
      CHECK(matrix_2.get_xy() == doctest::Approx(5.75));
      CHECK(matrix_2.get_xz() == doctest::Approx(-8.25));
      CHECK(matrix_2.get_yx() == doctest::Approx(10.5));
      CHECK(matrix_2.get_yy() == doctest::Approx(-1.25));
      CHECK(matrix_2.get_yz() == doctest::Approx(4));
      CHECK(matrix_2.get_zx() == doctest::Approx(-4.25));
      CHECK(matrix_2.get_zy() == doctest::Approx(9.5));
      CHECK(matrix_2.get_zz() == doctest::Approx(4.75));
    }
  }
  SUBCASE("test get/set functions") {
    SUBCASE("test get functions") {
      math::Matrix3x3 matrix_2 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      CHECK(matrix_2.get_xx() == doctest::Approx(54));
      CHECK(matrix_2.get_xy() == doctest::Approx(23));
      CHECK(matrix_2.get_xz() == doctest::Approx(-33));
      CHECK(matrix_2.get_yx() == doctest::Approx(42));
      CHECK(matrix_2.get_yy() == doctest::Approx(-5));
      CHECK(matrix_2.get_yz() == doctest::Approx(16));
      CHECK(matrix_2.get_zx() == doctest::Approx(-17));
      CHECK(matrix_2.get_zy() == doctest::Approx(38));
      CHECK(matrix_2.get_zz() == doctest::Approx(19));
    }
    SUBCASE("test set functions") {
      math::Matrix3x3 matrix_1 = {};
      matrix_1.set_xx(1);
      matrix_1.set_xy(1);
      matrix_1.set_xz(1);
      matrix_1.set_yx(1);
      matrix_1.set_yy(1);
      matrix_1.set_yz(1);
      matrix_1.set_zx(1);
      matrix_1.set_zy(1);
      matrix_1.set_zz(1);
      CHECK(matrix_1.get_xx() == doctest::Approx(1));
      CHECK(matrix_1.get_xy() == doctest::Approx(1));
      CHECK(matrix_1.get_xz() == doctest::Approx(1));
      CHECK(matrix_1.get_yx() == doctest::Approx(1));
      CHECK(matrix_1.get_yy() == doctest::Approx(1));
      CHECK(matrix_1.get_yz() == doctest::Approx(1));
      CHECK(matrix_1.get_zx() == doctest::Approx(1));
      CHECK(matrix_1.get_zy() == doctest::Approx(1));
      CHECK(matrix_1.get_zz() == doctest::Approx(1));
    }
  }
  SUBCASE("test matrix/matrix math") {
    SUBCASE("test matrix/matrix addition") {
      SUBCASE("test zero addition") {
        SUBCASE("test 0/0 addition") {
          math::Matrix3x3 matrix_1 = {};
          math::Matrix3x3 matrix_2 = {};
          math::Matrix3x3 matrix_3 = matrix_1 + matrix_2;
          CHECK(matrix_3.get_xx() == doctest::Approx(0));
          CHECK(matrix_3.get_xy() == doctest::Approx(0));
          CHECK(matrix_3.get_xz() == doctest::Approx(0));
          CHECK(matrix_3.get_yx() == doctest::Approx(0));
          CHECK(matrix_3.get_yy() == doctest::Approx(0));
          CHECK(matrix_3.get_yz() == doctest::Approx(0));
          CHECK(matrix_3.get_zx() == doctest::Approx(0));
          CHECK(matrix_3.get_zy() == doctest::Approx(0));
          CHECK(matrix_3.get_zz() == doctest::Approx(0));
        }
        SUBCASE("test 0/x addition") {
          math::Matrix3x3 matrix_1 = {};
          math::Matrix3x3 matrix_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
          math::Matrix3x3 matrix_3 = matrix_1 + matrix_2;
          CHECK(matrix_3.get_xx() == doctest::Approx(0));
          CHECK(matrix_3.get_xy() == doctest::Approx(1));
          CHECK(matrix_3.get_xz() == doctest::Approx(2));
          CHECK(matrix_3.get_yx() == doctest::Approx(3));
          CHECK(matrix_3.get_yy() == doctest::Approx(4));
          CHECK(matrix_3.get_yz() == doctest::Approx(5));
          CHECK(matrix_3.get_zx() == doctest::Approx(6));
          CHECK(matrix_3.get_zy() == doctest::Approx(7));
          CHECK(matrix_3.get_zz() == doctest::Approx(8));
        }
      }
      SUBCASE("test simple addition") {
        math::Matrix3x3 matrix_1 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_3 = matrix_2 + matrix_1;
        CHECK(matrix_3.get_xx() == doctest::Approx(0));
        CHECK(matrix_3.get_xy() == doctest::Approx(2));
        CHECK(matrix_3.get_xz() == doctest::Approx(4));
        CHECK(matrix_3.get_yx() == doctest::Approx(6));
        CHECK(matrix_3.get_yy() == doctest::Approx(8));
        CHECK(matrix_3.get_yz() == doctest::Approx(10));
        CHECK(matrix_3.get_zx() == doctest::Approx(12));
        CHECK(matrix_3.get_zy() == doctest::Approx(14));
        CHECK(matrix_3.get_zz() == doctest::Approx(16));
      }
      SUBCASE("test complex addition") {
        math::Matrix3x3 matrix_1 = {12, 3.5, -8, 0, -3, 0.7, 22, 11, 90};
        math::Matrix3x3 matrix_2 = {-2, 23.5, 18, -1, 9, 4.1213, 2.71828, -1.41, 0};
        math::Matrix3x3 matrix_3 = matrix_1 + matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(10));
        CHECK(matrix_3.get_xy() == doctest::Approx(27));
        CHECK(matrix_3.get_xz() == doctest::Approx(10));
        CHECK(matrix_3.get_yx() == doctest::Approx(-1));
        CHECK(matrix_3.get_yy() == doctest::Approx(6));
        CHECK(matrix_3.get_yz() == doctest::Approx(4.8213));
        CHECK(matrix_3.get_zx() == doctest::Approx(24.71828));
        CHECK(matrix_3.get_zy() == doctest::Approx(9.59));
        CHECK(matrix_3.get_zz() == doctest::Approx(90));
      }
    }
    SUBCASE("test matrix subtraction") {
      SUBCASE("test zero subtraction") {
        math::Matrix3x3 matrix_1 = {1, 0, 1, -1, -1, 0, 1, -1, 0};
        math::Matrix3x3 matrix_2 = {1, 1, 0, -1, 0, -1, -1, 1, 0};
        math::Matrix3x3 matrix_3 = matrix_1 - matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(0));
        CHECK(matrix_3.get_xy() == doctest::Approx(-1));
        CHECK(matrix_3.get_xz() == doctest::Approx(1));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(-1));
        CHECK(matrix_3.get_yz() == doctest::Approx(1));
        CHECK(matrix_3.get_zx() == doctest::Approx(2));
        CHECK(matrix_3.get_zy() == doctest::Approx(-2));
        CHECK(matrix_3.get_zz() == doctest::Approx(0));
      }
      SUBCASE("test simple subtraction") {
        math::Matrix3x3 matrix_1 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_3 = matrix_1 - matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(0));
        CHECK(matrix_3.get_xy() == doctest::Approx(0));
        CHECK(matrix_3.get_xz() == doctest::Approx(0));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(0));
        CHECK(matrix_3.get_yz() == doctest::Approx(0));
        CHECK(matrix_3.get_zx() == doctest::Approx(0));
        CHECK(matrix_3.get_zy() == doctest::Approx(0));
        CHECK(matrix_3.get_zz() == doctest::Approx(0));
      }
      SUBCASE("test complex subtraction") {
        math::Matrix3x3 matrix_1 = {12, 3.5, -8, 0, -3, 0.7, 22, 11, 90};
        math::Matrix3x3 matrix_2 = {-2, 23.5, 18, -1, 9, 4.1213, 2.71828, -1.41, 0};
        math::Matrix3x3 matrix_3 = matrix_1 - matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(14));
        CHECK(matrix_3.get_xy() == doctest::Approx(-20));
        CHECK(matrix_3.get_xz() == doctest::Approx(-26));
        CHECK(matrix_3.get_yx() == doctest::Approx(1));
        CHECK(matrix_3.get_yy() == doctest::Approx(-12));
        CHECK(matrix_3.get_yz() == doctest::Approx(-3.4213));
        CHECK(matrix_3.get_zx() == doctest::Approx(19.28172));
        CHECK(matrix_3.get_zy() == doctest::Approx(12.41));
        CHECK(matrix_3.get_zz() == doctest::Approx(90));
      }
    }
    SUBCASE("test matrix multiplication") {
      SUBCASE("test zero multiplication") {
        math::Matrix3x3 matrix_1 = {};
        math::Matrix3x3 matrix_2 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_3 = matrix_1 * matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(0));
        CHECK(matrix_3.get_xy() == doctest::Approx(0));
        CHECK(matrix_3.get_xz() == doctest::Approx(0));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(0));
        CHECK(matrix_3.get_yz() == doctest::Approx(0));
        CHECK(matrix_3.get_zx() == doctest::Approx(0));
        CHECK(matrix_3.get_zy() == doctest::Approx(0));
        CHECK(matrix_3.get_zz() == doctest::Approx(0));
      }
      SUBCASE("test simple multiplication") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_3 = matrix_1 * matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(30));
        CHECK(matrix_3.get_xy() == doctest::Approx(36));
        CHECK(matrix_3.get_xz() == doctest::Approx(42));
        CHECK(matrix_3.get_yx() == doctest::Approx(66));
        CHECK(matrix_3.get_yy() == doctest::Approx(81));
        CHECK(matrix_3.get_yz() == doctest::Approx(96));
        CHECK(matrix_3.get_zx() == doctest::Approx(102));
        CHECK(matrix_3.get_zy() == doctest::Approx(126));
        CHECK(matrix_3.get_zz() == doctest::Approx(150));
      }
      SUBCASE("test complex multiplication") {
        math::Matrix3x3 matrix_1 = {12, 3.5, -8, 0, -3, 0.7, 22, 11, 90};
        math::Matrix3x3 matrix_2 = {-2, 23.5, 18, -1, 9, 4.1213, 2.71828, -1.41, 0};
        math::Matrix3x3 matrix_3 = matrix_1 * matrix_2;
        CHECK(matrix_3.get_xx() == doctest::Approx(-49.2462));
        CHECK(matrix_3.get_xy() == doctest::Approx(324.78));
        CHECK(matrix_3.get_xz() == doctest::Approx(230.425));
        CHECK(matrix_3.get_yx() == doctest::Approx(4.9028));
        CHECK(matrix_3.get_yy() == doctest::Approx(-27.987));
        CHECK(matrix_3.get_yz() == doctest::Approx(-12.3639));
        CHECK(matrix_3.get_zx() == doctest::Approx(189.645));
        CHECK(matrix_3.get_zy() == doctest::Approx(489.1));
        CHECK(matrix_3.get_zz() == doctest::Approx(441.334));
      }
      SUBCASE("test identity matrix multiplication") {
        math::Matrix3x3 matrix_1 = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        math::Matrix3x3 matrix_3 = {};
        math::Matrix3x3 matrix_2 = matrix_3.identity();
        math::Matrix3x3 matrix_4 = matrix_2 * matrix_1;
        CHECK(matrix_4.get_xx() == doctest::Approx(0));
        CHECK(matrix_4.get_xy() == doctest::Approx(1));
        CHECK(matrix_4.get_xz() == doctest::Approx(2));
        CHECK(matrix_4.get_yx() == doctest::Approx(3));
        CHECK(matrix_4.get_yy() == doctest::Approx(4));
        CHECK(matrix_4.get_yz() == doctest::Approx(5));
        CHECK(matrix_4.get_zx() == doctest::Approx(6));
        CHECK(matrix_4.get_zy() == doctest::Approx(7));
        CHECK(matrix_4.get_zz() == doctest::Approx(8));
      }
      SUBCASE("test zero multiplication void fxn") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {};
        math::Matrix3x3 matrix_3 = {};
        matrix_1.dot(matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(0));
        CHECK(matrix_3.get_xy() == doctest::Approx(0));
        CHECK(matrix_3.get_xz() == doctest::Approx(0));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(0));
        CHECK(matrix_3.get_yz() == doctest::Approx(0));
        CHECK(matrix_3.get_zx() == doctest::Approx(0));
        CHECK(matrix_3.get_zy() == doctest::Approx(0));
        CHECK(matrix_3.get_zz() == doctest::Approx(0));
      }
      SUBCASE("test simple multiplication void fxn") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_3 = {};
        matrix_1.dot(matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(30));
        CHECK(matrix_3.get_xy() == doctest::Approx(36));
        CHECK(matrix_3.get_xz() == doctest::Approx(42));
        CHECK(matrix_3.get_yx() == doctest::Approx(66));
        CHECK(matrix_3.get_yy() == doctest::Approx(81));
        CHECK(matrix_3.get_yz() == doctest::Approx(96));
        CHECK(matrix_3.get_zx() == doctest::Approx(102));
        CHECK(matrix_3.get_zy() == doctest::Approx(126));
        CHECK(matrix_3.get_zz() == doctest::Approx(150));
      }
      SUBCASE("test complex multiplication void fxn") {
        math::Matrix3x3 matrix_1 = {12, 3.5, -8, 0, -3, 0.7, 22, 11, 90};
        math::Matrix3x3 matrix_2 = {-2, 23.5, 18, -1, 9, 4.1213, 2.71828, -1.41, 0};
        math::Matrix3x3 matrix_3 = {};
        matrix_1.dot(matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(-49.2462));
        CHECK(matrix_3.get_xy() == doctest::Approx(324.78));
        CHECK(matrix_3.get_xz() == doctest::Approx(230.425));
        CHECK(matrix_3.get_yx() == doctest::Approx(4.9028));
        CHECK(matrix_3.get_yy() == doctest::Approx(-27.987));
        CHECK(matrix_3.get_yz() == doctest::Approx(-12.3639));
        CHECK(matrix_3.get_zx() == doctest::Approx(189.645));
        CHECK(matrix_3.get_zy() == doctest::Approx(489.1));
        CHECK(matrix_3.get_zz() == doctest::Approx(441.334));
      }
      SUBCASE("test identity multiplication void fxn") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = {};
        matrix_1.dot(matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(1));
        CHECK(matrix_3.get_xy() == doctest::Approx(2));
        CHECK(matrix_3.get_xz() == doctest::Approx(3));
        CHECK(matrix_3.get_yx() == doctest::Approx(4));
        CHECK(matrix_3.get_yy() == doctest::Approx(5));
        CHECK(matrix_3.get_yz() == doctest::Approx(6));
        CHECK(matrix_3.get_zx() == doctest::Approx(7));
        CHECK(matrix_3.get_zy() == doctest::Approx(8));
        CHECK(matrix_3.get_zz() == doctest::Approx(9));
      }
    }
    SUBCASE("test for presence of identity matrix") {
      math::Matrix3x3 matrix_1 = {};
      math::Matrix3x3 matrix_2 = matrix_1.identity();
      CHECK(matrix_2.get_xx() == doctest::Approx(1));
      CHECK(matrix_2.get_xy() == doctest::Approx(0));
      CHECK(matrix_2.get_xz() == doctest::Approx(0));
      CHECK(matrix_2.get_yx() == doctest::Approx(0));
      CHECK(matrix_2.get_yy() == doctest::Approx(1));
      CHECK(matrix_2.get_yz() == doctest::Approx(0));
      CHECK(matrix_2.get_zx() == doctest::Approx(0));
      CHECK(matrix_2.get_zy() == doctest::Approx(0));
      CHECK(matrix_2.get_zz() == doctest::Approx(1));
    }
    SUBCASE("test transposition functions") {
      SUBCASE("test transposition") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = matrix_1.transpose();
        CHECK(matrix_2.get_xx() == doctest::Approx(1));
        CHECK(matrix_2.get_xy() == doctest::Approx(4));
        CHECK(matrix_2.get_xz() == doctest::Approx(7));
        CHECK(matrix_2.get_yx() == doctest::Approx(2));
        CHECK(matrix_2.get_yy() == doctest::Approx(5));
        CHECK(matrix_2.get_yz() == doctest::Approx(8));
        CHECK(matrix_2.get_zx() == doctest::Approx(3));
        CHECK(matrix_2.get_zy() == doctest::Approx(6));
        CHECK(matrix_2.get_zz() == doctest::Approx(9));
      }
      SUBCASE("test transpoisition into matrix") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {};
        matrix_2.transpose(matrix_1);
        CHECK(matrix_2.get_xx() == doctest::Approx(1));
        CHECK(matrix_2.get_xy() == doctest::Approx(4));
        CHECK(matrix_2.get_xz() == doctest::Approx(7));
        CHECK(matrix_2.get_yx() == doctest::Approx(2));
        CHECK(matrix_2.get_yy() == doctest::Approx(5));
        CHECK(matrix_2.get_yz() == doctest::Approx(8));
        CHECK(matrix_2.get_zx() == doctest::Approx(3));
        CHECK(matrix_2.get_zy() == doctest::Approx(6));
        CHECK(matrix_2.get_zz() == doctest::Approx(9));
      }
      SUBCASE ("test get transposed") {
        math::Matrix3x3 matrix_0 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_1 = matrix_0.get_transposed();
        CHECK(matrix_1.get_xx() == doctest::Approx(1));
        CHECK(matrix_1.get_xy() == doctest::Approx(4));
        CHECK(matrix_1.get_xz() == doctest::Approx(7));
        CHECK(matrix_1.get_yx() == doctest::Approx(2));
        CHECK(matrix_1.get_yy() == doctest::Approx(5));
        CHECK(matrix_1.get_yz() == doctest::Approx(8));
        CHECK(matrix_1.get_zx() == doctest::Approx(3));
        CHECK(matrix_1.get_zy() == doctest::Approx(6));
        CHECK(matrix_1.get_zz() == doctest::Approx(9));
      }
    }
    SUBCASE("test matrix/vector multiplication") {
      SUBCASE("test zero multiplication") {
        math::Vector3 vec_1 = {0, 0, 0};
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Vector3 dotted_matrix_vector = matrix_1.dot(vec_1);// multiply_matrices(matrix_1, vec_1);
        CHECK(dotted_matrix_vector.get_x() == doctest::Approx(0));
        CHECK(dotted_matrix_vector.get_y() == doctest::Approx(0));
        CHECK(dotted_matrix_vector.get_z() == doctest::Approx(0));
      }
      SUBCASE("test simple multiplication") {
        math::Vector3 vec_1 = {1, 2, 3};
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Vector3 dotted_matrix_vector = matrix_1.dot(vec_1);
        CHECK(dotted_matrix_vector.get_x() == doctest::Approx(14));
        CHECK(dotted_matrix_vector.get_y() == doctest::Approx(32));
        CHECK(dotted_matrix_vector.get_z() == doctest::Approx(50));
      }
      SUBCASE("test complex multiplication") {
        math::Vector3 vec_1 = {4, 2.5, -2};
        math::Matrix3x3 matrix_1 = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
        math::Vector3 dotted_matrix_vector = matrix_1.dot(vec_1);
        CHECK(dotted_matrix_vector.get_x() == doctest::Approx(149));
        CHECK(dotted_matrix_vector.get_y() == doctest::Approx(-47.586));
        CHECK(dotted_matrix_vector.get_z() == doctest::Approx(21));
      }
      SUBCASE("test matrix/vector void multiplication") {
        SUBCASE("test zero multiplication") {
          math::Vector3 vec_1 = {};
          math::Vector3 vec_2 = {};
          math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
          matrix_1.dot(vec_1, vec_2);
          CHECK(vec_2.get_x() == doctest::Approx(0));
          CHECK(vec_2.get_y() == doctest::Approx(0));
          CHECK(vec_2.get_z() == doctest::Approx(0));
        }
        SUBCASE("test simple multiplication") {
          math::Vector3 vec_1 = {1, 2, 3};
          math::Vector3 vec_2 = {};
          math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
          matrix_1.dot(vec_1, vec_2);
          CHECK(vec_2.get_x() == doctest::Approx(14));
          CHECK(vec_2.get_y() == doctest::Approx(32));
          CHECK(vec_2.get_z() == doctest::Approx(50));
        }
        SUBCASE("test complex multiplication") {
          math::Vector3 vec_1 = {4, 2.5, -2};
          math::Vector3 vec_2 = {};
          math::Matrix3x3 matrix_1 = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
          matrix_1.dot(vec_1, vec_2);
          CHECK(vec_2.get_x() == doctest::Approx(149));
          CHECK(vec_2.get_y() == doctest::Approx(-47.586));
          CHECK(vec_2.get_z() == doctest::Approx(21));
        }
      }
    }
    SUBCASE("test flip orientation") {
      SUBCASE("test simple flip orientation") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = matrix_1.get_flip_orientation();
        CHECK(matrix_2.get_xx() == doctest::Approx(1));
        CHECK(matrix_2.get_xy() == doctest::Approx(2));
        CHECK(matrix_2.get_xz() == doctest::Approx(3));
        CHECK(matrix_2.get_yx() == doctest::Approx(-4));
        CHECK(matrix_2.get_yy() == doctest::Approx(-5));
        CHECK(matrix_2.get_yz() == doctest::Approx(-6));
        CHECK(matrix_2.get_zx() == doctest::Approx(-7));
        CHECK(matrix_2.get_zy() == doctest::Approx(-8));
        CHECK(matrix_2.get_zz() == doctest::Approx(-9));
      }
      SUBCASE("test complex flip orientation") {
        math::Matrix3x3 matrix_1 = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
        math::Matrix3x3 matrix_2 = matrix_1.get_flip_orientation();
        CHECK(matrix_2.get_xx() == doctest::Approx(5));
        CHECK(matrix_2.get_xy() == doctest::Approx(18));
        CHECK(matrix_2.get_xz() == doctest::Approx(-42));
        CHECK(matrix_2.get_yx() == doctest::Approx(0));
        CHECK(matrix_2.get_yy() == doctest::Approx(19));
        CHECK(matrix_2.get_yz() == doctest::Approx(-0.043));
        CHECK(matrix_2.get_zx() == doctest::Approx(-32));
        CHECK(matrix_2.get_zy() == doctest::Approx(29.2));
        CHECK(matrix_2.get_zz() == doctest::Approx(-17));
      }
    }
    SUBCASE("test difference between matrices") {
      SUBCASE("test simple difference") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {2, 3, 4, 5, 6, 7, 8, 9, 0};
        CHECK(matrix_2.difference(matrix_1) == doctest::Approx(17));
      }
      SUBCASE("test complex difference") {
        math::Matrix3x3 matrix_1 = {12, 3.5, -8, 0, -3, 0.7, 22, 11, 90};
        math::Matrix3x3 matrix_2 = {-2, 23.5, 18, -1, 9, 4.1213, 2.71828, -1.41, 0};
        CHECK(matrix_2.difference(matrix_1) == doctest::Approx(198.11302));
      }
    }
  }
  SUBCASE("test stringifications") {
    SUBCASE("string to matrix") {
      SUBCASE("test normal circumtances") {
        String matrix_1_string = "1 2 3 4 5 6 7 8 9";
        math::Matrix3x3 matrix_1 = math::matrix_from_str(matrix_1_string);
        CHECK(matrix_1.get_xx() == doctest::Approx(1));
        CHECK(matrix_1.get_xy() == doctest::Approx(2));
        CHECK(matrix_1.get_xz() == doctest::Approx(3));
        CHECK(matrix_1.get_yx() == doctest::Approx(4));
        CHECK(matrix_1.get_yy() == doctest::Approx(5));
        CHECK(matrix_1.get_yz() == doctest::Approx(6));
        CHECK(matrix_1.get_zx() == doctest::Approx(7));
        CHECK(matrix_1.get_zy() == doctest::Approx(8));
        CHECK(matrix_1.get_zz() == doctest::Approx(9));
      }
      SUBCASE("test too many inputs") {
        String matrix_1_string = "1 2 3 4 5 6 7 8 9 0";
        CHECK_THROWS_AS(math::matrix_from_str(matrix_1_string), base::InputException);
      }
      SUBCASE("test too few inputs") {
        String matrix_1_string = "1 2 3 4 5 6 7 8";
        CHECK_THROWS_AS(math::matrix_from_str(matrix_1_string), base::InputException);
      }
      SUBCASE("test no inputs") {
        String matrix_1_string = "";
        CHECK_THROWS_AS(math::matrix_from_str(matrix_1_string), base::InputException);
      }
    }
  }
  SUBCASE("test transform 1") {
    SUBCASE("test simple transform") {
      math::Matrix3x3 matrix = {0, 1, 2, 3, 4, 5, 6, 7, 8};
      math::Matrix3x3 matrix_1 = matrix.transform_1();
      CHECK(matrix_1.get_xx() == doctest::Approx(0));
      CHECK(matrix_1.get_xy() == doctest::Approx(1));
      CHECK(matrix_1.get_xz() == doctest::Approx(2));
      CHECK(matrix_1.get_yx() == doctest::Approx(-3));
      CHECK(matrix_1.get_yy() == doctest::Approx(-4));
      CHECK(matrix_1.get_yz() == doctest::Approx(-5));
      CHECK(matrix_1.get_zx() == doctest::Approx(-6));
      CHECK(matrix_1.get_zy() == doctest::Approx(-7));
      CHECK(matrix_1.get_zz() == doctest::Approx(-8));
    }
    SUBCASE("test complex transform") {
      math::Matrix3x3 matrix = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
      math::Matrix3x3 matrix_1 = matrix.transform_1();
      CHECK(matrix_1.get_xx() == doctest::Approx(5));
      CHECK(matrix_1.get_xy() == doctest::Approx(18));
      CHECK(matrix_1.get_xz() == doctest::Approx(-42));
      CHECK(matrix_1.get_yx() == doctest::Approx(0));
      CHECK(matrix_1.get_yy() == doctest::Approx(19));
      CHECK(matrix_1.get_yz() == doctest::Approx(-0.043));
      CHECK(matrix_1.get_zx() == doctest::Approx(-32));
      CHECK(matrix_1.get_zy() == doctest::Approx(29.2));
      CHECK(matrix_1.get_zz() == doctest::Approx(-17));
    }
  }
  SUBCASE("test unitarize") {
    SUBCASE("test simple unitarize") {
      math::Matrix3x3 matrix_1 = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
      math::Matrix3x3 matrix_2 = matrix_1.get_unitarize();
      CHECK(matrix_2.get_xx() == doctest::Approx(float(1 / sqrt(14))));
      CHECK(matrix_2.get_xy() == doctest::Approx(float(sqrt(0.285714))));
      CHECK(matrix_2.get_xz() == doctest::Approx(float(3 / sqrt(14))));
      CHECK(matrix_2.get_yx() == doctest::Approx(float(4 / sqrt(77))));
      CHECK(matrix_2.get_yy() == doctest::Approx(float(5 / sqrt(77))));
      CHECK(matrix_2.get_yz() == doctest::Approx(float(6 / sqrt(77))));
      CHECK(matrix_2.get_zx() == doctest::Approx(float(7 / sqrt(194))));
      CHECK(matrix_2.get_zy() == doctest::Approx(0.5743665269));
      CHECK(matrix_2.get_zz() == doctest::Approx(float(9 / sqrt(194))));
      CHECK(sqrt((matrix_2.get_xx() * matrix_2.get_xx()) + (matrix_2.get_xy() * matrix_2.get_xy()) + (matrix_2.get_xz() * matrix_2.get_xz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_yx() * matrix_2.get_yx()) + (matrix_2.get_yy() * matrix_2.get_yy()) + (matrix_2.get_yz() * matrix_2.get_yz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_zx() * matrix_2.get_zx()) + (matrix_2.get_zy() * matrix_2.get_zy()) + (matrix_2.get_zz() * matrix_2.get_zz())) == doctest::Approx(1.0f));
    }
    SUBCASE("test complex unitarize") {
      math::Matrix3x3 matrix_1 = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
      math::Matrix3x3 matrix_2 = matrix_1.get_unitarize();
      CHECK(matrix_2.get_xx() == doctest::Approx(0.1087727869));
      CHECK(matrix_2.get_xy() == doctest::Approx(float(18/sqrt(2113))));
      CHECK(matrix_2.get_xz() == doctest::Approx(float(-42/sqrt(2113))));
      CHECK(matrix_2.get_yx() == doctest::Approx(0));
      CHECK(matrix_2.get_yy() == doctest::Approx(-0.999997));
      CHECK(matrix_2.get_yz() == doctest::Approx(2.26315e-3));
      CHECK(matrix_2.get_zx() == doctest::Approx(0.687633));
      CHECK(matrix_2.get_zy() == doctest::Approx(-0.627465));
      CHECK(matrix_2.get_zz() == doctest::Approx(0.365305));
      CHECK(sqrt((matrix_2.get_xx() * matrix_2.get_xx()) + (matrix_2.get_xy() * matrix_2.get_xy()) + (matrix_2.get_xz() * matrix_2.get_xz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_yx() * matrix_2.get_yx()) + (matrix_2.get_yy() * matrix_2.get_yy()) + (matrix_2.get_yz() * matrix_2.get_yz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_zx() * matrix_2.get_zx()) + (matrix_2.get_zy() * matrix_2.get_zy()) + (matrix_2.get_zz() * matrix_2.get_zz())) == doctest::Approx(1.0f));

    }
    SUBCASE("test void simple unitarize") {
      math::Matrix3x3 matrix_2 = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f};
      matrix_2.unitarize();
      CHECK(matrix_2.get_xx() == doctest::Approx(float(1 / sqrt(14))));
      CHECK(matrix_2.get_xy() == doctest::Approx(float(sqrt(0.285714))));
      CHECK(matrix_2.get_xz() == doctest::Approx(float(3 / sqrt(14))));
      CHECK(matrix_2.get_yx() == doctest::Approx(float(4 / sqrt(77))));
      CHECK(matrix_2.get_yy() == doctest::Approx(float(5 / sqrt(77))));
      CHECK(matrix_2.get_yz() == doctest::Approx(float(6 / sqrt(77))));
      CHECK(matrix_2.get_zx() == doctest::Approx(float(7 / sqrt(194))));
      CHECK(matrix_2.get_zy() == doctest::Approx(0.5743665269));
      CHECK(matrix_2.get_zz() == doctest::Approx(float(9 / sqrt(194))));
      CHECK(sqrt((matrix_2.get_xx() * matrix_2.get_xx()) + (matrix_2.get_xy() * matrix_2.get_xy()) + (matrix_2.get_xz() * matrix_2.get_xz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_yx() * matrix_2.get_yx()) + (matrix_2.get_yy() * matrix_2.get_yy()) + (matrix_2.get_yz() * matrix_2.get_yz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_zx() * matrix_2.get_zx()) + (matrix_2.get_zy() * matrix_2.get_zy()) + (matrix_2.get_zz() * matrix_2.get_zz())) == doctest::Approx(1.0f));
    }
    SUBCASE("test void complex unitarize") {
      math::Matrix3x3 matrix_2 = {5, 18, -42, 0, -19, .043, 32, -29.2, 17};
      matrix_2.unitarize();
      CHECK(matrix_2.get_xx() == doctest::Approx(0.1087727869));
      CHECK(matrix_2.get_xy() == doctest::Approx(float(18/sqrt(2113))));
      CHECK(matrix_2.get_xz() == doctest::Approx(float(-42/sqrt(2113))));
      CHECK(matrix_2.get_yx() == doctest::Approx(0));
      CHECK(matrix_2.get_yy() == doctest::Approx(-0.999997));
      CHECK(matrix_2.get_yz() == doctest::Approx(2.26315e-3));
      CHECK(matrix_2.get_zx() == doctest::Approx(0.687633));
      CHECK(matrix_2.get_zy() == doctest::Approx(-0.627465));
      CHECK(matrix_2.get_zz() == doctest::Approx(0.365305));
      CHECK(sqrt((matrix_2.get_xx() * matrix_2.get_xx()) + (matrix_2.get_xy() * matrix_2.get_xy()) + (matrix_2.get_xz() * matrix_2.get_xz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_yx() * matrix_2.get_yx()) + (matrix_2.get_yy() * matrix_2.get_yy()) + (matrix_2.get_yz() * matrix_2.get_yz())) == doctest::Approx(1.0f));
      CHECK(sqrt((matrix_2.get_zx() * matrix_2.get_zx()) + (matrix_2.get_zy() * matrix_2.get_zy()) + (matrix_2.get_zz() * matrix_2.get_zz())) == doctest::Approx(1.0f));

    }
  }
}

