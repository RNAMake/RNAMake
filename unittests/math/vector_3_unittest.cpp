//
// Created by Joe Yesselman on 4/18/22.
//

#include "../common.hpp"

#include <math/vector_3.hpp>

TEST_CASE("Test xyz vector ") {
  SUBCASE("test construction") {
    SUBCASE("test no args") {
      math::Vector3 vec{};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(0));
      CHECK(vec.get_z() == doctest::Approx(0));
    }
    SUBCASE("test supply ints") {
      math::Vector3 vec = {0, 1, 2};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(1));
      CHECK(vec.get_z() == doctest::Approx(2));
    }
    SUBCASE("test supply floats") {
      math::Vector3 vec = {0.0f, 1.0f, 2.0f};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(1));
      CHECK(vec.get_z() == doctest::Approx(2));
    }
    SUBCASE("from double vector") {
      Reals nums = {0.0, 1.0, 2.0};
      math::Vector3 vec(nums);
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(1));
      CHECK(vec.get_z() == doctest::Approx(2));
      // TODO check to make sure this actually raises an assert?
      Reals nums2 = {0.0, 1.0, 2.0, 3.0};
      math::Vector3 vec2(nums2);
    }
  }
  SUBCASE("test copy") {
    SUBCASE("test trival copy") {
      math::Vector3 vec = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec2{vec};
      vec.set_x(1.1);
      vec.set_y(1.1);
      vec.set_z(1.1);
      // should still have original values
      CHECK(vec2.get_x() == doctest::Approx(0));
      CHECK(vec2.get_y() == doctest::Approx(1));
      CHECK(vec2.get_z() == doctest::Approx(2));
    }
    SUBCASE("using = operator") {
      math::Vector3 vec = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_2 = vec;
      vec.set_x(1.1);
      vec.set_y(1.1);
      vec.set_z(1.1);
      // should still have original values
      CHECK(vec_2.get_x() == doctest::Approx(0));
      CHECK(vec_2.get_y() == doctest::Approx(1));
      CHECK(vec_2.get_z() == doctest::Approx(2));
    }
  }
  SUBCASE("test operators") {
    SUBCASE("test +") {
      SUBCASE("vec + vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_3 = vec_1 + vec_2;
        CHECK(vec_3.get_x() == doctest::Approx(0));
        CHECK(vec_3.get_y() == doctest::Approx(2));
        CHECK(vec_3.get_z() == doctest::Approx(4));
      }
      SUBCASE("vec + const") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real c = 0.1f;
        vec_1 = vec_1 + c;
        CHECK(vec_1.get_x() == doctest::Approx(0.1));
        CHECK(vec_1.get_y() == doctest::Approx(1.1));
        CHECK(vec_1.get_z() == doctest::Approx(2.1));
      }
      SUBCASE("vec + const limit") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real c = 0.0001f;
        vec_1 = vec_1 + c;
        CHECK(vec_1.get_x() == doctest::Approx(0.0001));
        CHECK(vec_1.get_y() == doctest::Approx(1.0001));
        CHECK(vec_1.get_z() == doctest::Approx(2.0001));
      }
    }
    SUBCASE("test -") {
      SUBCASE("vec - vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_3 = vec_1 - vec_2;
        CHECK(vec_3.get_x() == doctest::Approx(0));
        CHECK(vec_3.get_y() == doctest::Approx(0));
        CHECK(vec_3.get_z() == doctest::Approx(0));
      }
      SUBCASE("vec - constant") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        vec_1 = vec_1 - 3;
        CHECK(vec_1.get_x() == doctest::Approx(-3.0));
        CHECK(vec_1.get_y() == doctest::Approx(-2.0));
        CHECK(vec_1.get_z() == doctest::Approx(-1.0));
      }
    }
    SUBCASE("test *") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_2 = vec_1 * 2;
      CHECK(vec_2.get_x() == doctest::Approx(0));
      CHECK(vec_2.get_y() == doctest::Approx(2));
      CHECK(vec_2.get_z() == doctest::Approx(4));
    }
    SUBCASE("test /") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_2 = vec_1 / 2;
      CHECK(vec_2.get_x() == doctest::Approx(0));
      CHECK(vec_2.get_y() == doctest::Approx(0.5));
      CHECK(vec_2.get_z() == doctest::Approx(1));
    }
    SUBCASE("test / 0") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      CHECK_THROWS_AS(vec_1 / 0, base::MathException);
    }
    SUBCASE("test +=") {
      SUBCASE("test vec += vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        vec_1 += vec_2;
        CHECK(vec_1.get_x() == doctest::Approx(0));
        CHECK(vec_1.get_y() == doctest::Approx(2));
        CHECK(vec_1.get_z() == doctest::Approx(4));
      }
      SUBCASE("test vec += constant") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        vec_1 += 1.0f;
        CHECK(vec_1.get_x() == doctest::Approx(1));
        CHECK(vec_1.get_y() == doctest::Approx(2));
        CHECK(vec_1.get_z() == doctest::Approx(3));
      }
    }
  }
}
