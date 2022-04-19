//
// Created by Joe Yesselman on 4/18/22.
//

#include "../common.hpp"

#include <math/xyz_vector.hpp>

TEST_CASE("Test xyz vector ") {
  SUBCASE("test construction") {
    SUBCASE("test no args") {
      math::xyzVector vec{};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(0));
      CHECK(vec.get_z() == doctest::Approx(0));
    }
    SUBCASE("test supply ints") {
      math::xyzVector vec = {0, 1, 2};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(1));
      CHECK(vec.get_z() == doctest::Approx(2));
    }
    SUBCASE("test supply floats") {
      math::xyzVector vec = {0.0f, 1.0f, 2.0f};
      CHECK(vec.get_x() == doctest::Approx(0));
      CHECK(vec.get_y() == doctest::Approx(1));
      CHECK(vec.get_z() == doctest::Approx(2));
    }
  }
  SUBCASE("test copy") {
    SUBCASE("test trival copy") {
      math::xyzVector vec = {0.0f, 1.0f, 2.0f};
      math::xyzVector vec2{vec};
      vec.set_x(1.1);
      vec.set_y(1.1);
      vec.set_z(1.1);
      // should still have original values
      CHECK(vec2.get_x() == doctest::Approx(0));
      CHECK(vec2.get_y() == doctest::Approx(1));
      CHECK(vec2.get_z() == doctest::Approx(2));
    }
    SUBCASE("using = operator") {
      math::xyzVector vec = {0.0f, 1.0f, 2.0f};
      math::xyzVector vec2 = vec;
      vec.set_x(1.1);
      vec.set_y(1.1);
      vec.set_z(1.1);
      // should still have original values
      CHECK(vec2.get_x() == doctest::Approx(0));
      CHECK(vec2.get_y() == doctest::Approx(1));
      CHECK(vec2.get_z() == doctest::Approx(2));
    }
  }
}