//
// Created by Joe Yesselman on 6/8/22.
//

#include "../common.hpp"
#include <math/rotation.hpp>
#include <math.h>

TEST_CASE("test rotation math ") {
  SUBCASE("test rotation between frames") {
    SUBCASE("test simple rotation between two frames") {
      SUBCASE("test reflection around x-axis") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5,6 ,7 , 8, 9};
        math::Matrix3x3 matrix_2 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(-0.123091491));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.4923659639));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.8616404369));
        CHECK(matrix_3.get_yx() == doctest::Approx(-0.2073903389));
        CHECK(matrix_3.get_yy() == doctest::Approx(0.5184758474));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.8295613558));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.2672612419));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.5345224838));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.8017837257));
      }
      SUBCASE("test reflection around y-axis") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, -1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.123091491));
        CHECK(matrix_3.get_xy() == doctest::Approx(-0.4923659639));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.8616404369));
        CHECK(matrix_3.get_yx() == doctest::Approx(0.2073903389));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.5184758474));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.8295613558));
        CHECK(matrix_3.get_zx() == doctest::Approx(0.2672612419));
        CHECK(matrix_3.get_zy() == doctest::Approx(-0.5345224838));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.8017837257));
      }
      SUBCASE("test reflection around z-axis") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.123091491));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.4923659639));
        CHECK(matrix_3.get_xz() == doctest::Approx(-0.8616404369));
        CHECK(matrix_3.get_yx() == doctest::Approx(0.2073903389));
        CHECK(matrix_3.get_yy() == doctest::Approx(0.5184758474));
        CHECK(matrix_3.get_yz() == doctest::Approx(-0.8295613558));
        CHECK(matrix_3.get_zx() == doctest::Approx(0.2672612419));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.5345224838));
        CHECK(matrix_3.get_zz() == doctest::Approx(-0.8017837257));
      }
    }
    SUBCASE("test complex rotation between two frames") {
      SUBCASE("test reflection around x-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(-0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(-0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.4600066262));
      }
      SUBCASE("test reflection around y-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, -1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(-0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(-0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(-0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.4600066262));
      }
      SUBCASE("test reflection around z-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        math::Matrix3x3 matrix_3 = math::rotation_between_frames(matrix_1, matrix_2);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(-0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(-0.4600066262));
      }
    }
  }
  SUBCASE("test void rotation between frames") {
    SUBCASE("test simple rotation between two frames") {
      SUBCASE("test reflection around x-axis") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6,7 , 8, 9};
        math::Matrix3x3 matrix_2 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = {};
        math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(-0.123091491));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.4923659639));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.8616404369));
        CHECK(matrix_3.get_yx() == doctest::Approx(-0.2073903389));
        CHECK(matrix_3.get_yy() == doctest::Approx(0.5184758474));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.8295613558));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.2672612419));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.5345224838));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.8017837257));
      }
      SUBCASE("test reflection around y-axis") {
        math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, -1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = {};
        math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.123091491));
        CHECK(matrix_3.get_xy() == doctest::Approx(-0.4923659639));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.8616404369));
        CHECK(matrix_3.get_yx() == doctest::Approx(0.2073903389));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.5184758474));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.8295613558));
        CHECK(matrix_3.get_zx() == doctest::Approx(0.2672612419));
        CHECK(matrix_3.get_zy() == doctest::Approx(-0.5345224838));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.8017837257));
      }
      SUBCASE("test reflection around z-axis") {
        SUBCASE("test reflection around z-axis") {
          math::Matrix3x3 matrix_1 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
          math::Matrix3x3 matrix_2 = {1, 0, 0, 0, 1, 0, 0, 0, -1};
          math::Matrix3x3 matrix_3 = {};
          math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
          CHECK(matrix_3.get_xx() == doctest::Approx(0.123091491));
          CHECK(matrix_3.get_xy() == doctest::Approx(0.4923659639));
          CHECK(matrix_3.get_xz() == doctest::Approx(-0.8616404369));
          CHECK(matrix_3.get_yx() == doctest::Approx(0.2073903389));
          CHECK(matrix_3.get_yy() == doctest::Approx(0.5184758474));
          CHECK(matrix_3.get_yz() == doctest::Approx(-0.8295613558));
          CHECK(matrix_3.get_zx() == doctest::Approx(0.2672612419));
          CHECK(matrix_3.get_zy() == doctest::Approx(0.5345224838));
          CHECK(matrix_3.get_zz() == doctest::Approx(-0.8017837257));
        }
      }
    }
    SUBCASE("test complex rotation between two frames") {
      SUBCASE("test reflection around x-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {-1, 0, 0, 0, 1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = {};
        math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(-0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(-0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.4600066262));
      }
      SUBCASE("test reflection around y-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, -1, 0, 0, 0, 1};
        math::Matrix3x3 matrix_3 = {};
        math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(-0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(-0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(-0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(0.4600066262));
      }
      SUBCASE("test reflection around z-axis") {
        math::Matrix3x3 matrix_1 = {54, 0, -33, 42, -5, 16, -17, 38, 19};
        math::Matrix3x3 matrix_2 = {1, 0, 0, 0, 1, 0, 0, 0, -1};
        math::Matrix3x3 matrix_3 = {};
        math::rotation_between_frames(matrix_1, matrix_2, matrix_3);
        CHECK(matrix_3.get_xx() == doctest::Approx(0.7660537828));
        CHECK(matrix_3.get_xy() == doctest::Approx(0.5958196088));
        CHECK(matrix_3.get_xz() == doctest::Approx(0.2411650798));
        CHECK(matrix_3.get_yx() == doctest::Approx(0));
        CHECK(matrix_3.get_yy() == doctest::Approx(-0.1304545126));
        CHECK(matrix_3.get_yz() == doctest::Approx(-0.9914542955));
        CHECK(matrix_3.get_zx() == doctest::Approx(-0.7989588771));
        CHECK(matrix_3.get_zy() == doctest::Approx(0.3873740010));
        CHECK(matrix_3.get_zz() == doctest::Approx(-0.4600066262));
      }
    }
  }
  SUBCASE("test difference between two rotation matrices") {
    // TODO be implemented
  }
}