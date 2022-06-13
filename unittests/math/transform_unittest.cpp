//
// Created by Joe Yesselman on 6/10/22.
//

#include "../common.hpp"
#include <sstream>
#include <math/transform.hpp>

TEST_CASE ("test rotation math ") {
  SUBCASE("test constructors") {
    SUBCASE("test empty constructor") {
      math::Transform transform = math::Transform();
      CHECK(transform.get_xx() == doctest::Approx(1));
      CHECK(transform.get_xy() == doctest::Approx(0));
      CHECK(transform.get_xz() == doctest::Approx(0));
      CHECK(transform.get_yx() == doctest::Approx(0));
      CHECK(transform.get_yy() == doctest::Approx(1));
      CHECK(transform.get_yz() == doctest::Approx(0));
      CHECK(transform.get_zx() == doctest::Approx(0));
      CHECK(transform.get_zy() == doctest::Approx(0));
      CHECK(transform.get_zz() == doctest::Approx(1));
      CHECK(transform.get_px() == doctest::Approx(0));
      CHECK(transform.get_py() == doctest::Approx(0));
      CHECK(transform.get_pz() == doctest::Approx(0));
    }
    SUBCASE("test zero constructor") {
      math::Transform transform = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      CHECK(transform.get_xx() == doctest::Approx(0));
      CHECK(transform.get_xy() == doctest::Approx(0));
      CHECK(transform.get_xz() == doctest::Approx(0));
      CHECK(transform.get_yx() == doctest::Approx(0));
      CHECK(transform.get_yy() == doctest::Approx(0));
      CHECK(transform.get_yz() == doctest::Approx(0));
      CHECK(transform.get_zx() == doctest::Approx(0));
      CHECK(transform.get_zy() == doctest::Approx(0));
      CHECK(transform.get_zz() == doctest::Approx(0));
      CHECK(transform.get_px() == doctest::Approx(0));
      CHECK(transform.get_py() == doctest::Approx(0));
      CHECK(transform.get_pz() == doctest::Approx(0));
    }
    SUBCASE("test number constructor") {
      math::Transform transform = {1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3};
      CHECK(transform.get_xx() == doctest::Approx(1));
      CHECK(transform.get_xy() == doctest::Approx(2));
      CHECK(transform.get_xz() == doctest::Approx(3));
      CHECK(transform.get_yx() == doctest::Approx(4));
      CHECK(transform.get_yy() == doctest::Approx(5));
      CHECK(transform.get_yz() == doctest::Approx(6));
      CHECK(transform.get_zx() == doctest::Approx(7));
      CHECK(transform.get_zy() == doctest::Approx(8));
      CHECK(transform.get_zz() == doctest::Approx(9));
      CHECK(transform.get_px() == doctest::Approx(1));
      CHECK(transform.get_py() == doctest::Approx(2));
      CHECK(transform.get_pz() == doctest::Approx(3));
    }
    SUBCASE("test complicated number constructor") {
      math::Transform transform = {5, -12, 0.001, 6, 0, -2, 2.71828, 1.41, 9, 1, 0, 0};
      CHECK(transform.get_xx() == doctest::Approx(5));
      CHECK(transform.get_xy() == doctest::Approx(-12));
      CHECK(transform.get_xz() == doctest::Approx(0.001));
      CHECK(transform.get_yx() == doctest::Approx(6));
      CHECK(transform.get_yy() == doctest::Approx(0));
      CHECK(transform.get_yz() == doctest::Approx(-2));
      CHECK(transform.get_zx() == doctest::Approx(2.71828));
      CHECK(transform.get_zy() == doctest::Approx(1.41));
      CHECK(transform.get_zz() == doctest::Approx(9));
      CHECK(transform.get_px() == doctest::Approx(1));
      CHECK(transform.get_py() == doctest::Approx(0));
      CHECK(transform.get_pz() == doctest::Approx(0));
    }
  }
  SUBCASE("test void dot function") {
    SUBCASE("test simple void dot function") {
      math::Transform matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0};
      math::Matrix3x3 matrix_2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
      math::Transform transform_1 = {};
      matrix_1.dot(matrix_2, transform_1);
      CHECK(transform_1.get_xx() == doctest::Approx(30));
      CHECK(transform_1.get_xy() == doctest::Approx(36));
      CHECK(transform_1.get_xz() == doctest::Approx(42));
      CHECK(transform_1.get_yx() == doctest::Approx(66));
      CHECK(transform_1.get_yy() == doctest::Approx(81));
      CHECK(transform_1.get_yz() == doctest::Approx(96));
      CHECK(transform_1.get_zx() == doctest::Approx(102));
      CHECK(transform_1.get_zy() == doctest::Approx(126));
      CHECK(transform_1.get_zz() == doctest::Approx(150));
      CHECK(transform_1.get_px() == doctest::Approx(0));
      CHECK(transform_1.get_py() == doctest::Approx(0));
      CHECK(transform_1.get_pz() == doctest::Approx(0));
    }
    SUBCASE("test complex void dot function") {
      math::Transform transform_1 = {5, -12, 0.001, 6, 0, -2, 2.71828, 1.41, 9, 1.0f, 0, 6.28f};
      math::Matrix3x3 matrix_1 = {54, 23, -33, 42, -5, 16, -17, 38, 19};
      math::Transform transform_2 = {};
      transform_1.dot(matrix_1, transform_2);
      CHECK(transform_2.get_xx() == doctest::Approx(-234.017));
      CHECK(transform_2.get_xy() == doctest::Approx(175.038));
      CHECK(transform_2.get_xz() == doctest::Approx(-356.981));
      CHECK(transform_2.get_yx() == doctest::Approx(358));
      CHECK(transform_2.get_yy() == doctest::Approx(62));
      CHECK(transform_2.get_yz() == doctest::Approx(-236));
      CHECK(transform_2.get_zx() == doctest::Approx(53.00712));
      CHECK(transform_2.get_zy() == doctest::Approx(397.47044));
      CHECK(transform_2.get_zz() == doctest::Approx(103.85676));
      // TODO fix this
      //CHECK(transform_2.get_px() == doctest::Approx(0));
      //CHECK(transform_2.get_py() == doctest::Approx(0));
      //CHECK(transform_2.get_pz() == doctest::Approx(0));
    }
  }
}