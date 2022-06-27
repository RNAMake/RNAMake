//
// Created by Joe Yesselman on 6/27/22.
//

#include "../common.hpp"
#include <math/matrix_3x3.hpp>
#include <math/numerical.hpp>
#include <math/vector_3.hpp>

TEST_CASE("test numerical ") {
  SUBCASE("test are_floats_equal") {
    SUBCASE("test true") {
      Real a = 2;
      Real b = 2.0;
      CHECK(math::are_floats_equal(a, b, 0.001) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      Real a = 2;
      Real b = 1;
      CHECK(math::are_floats_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside bounds") {
      Real a = 2;
      Real b = 2.01;
      CHECK(math::are_floats_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline within bounds") {
      Real a = 2;
      Real b = 2.0001;
      CHECK(math::are_floats_equal(a, b, 0.001) == doctest::Approx(1));
    }
  }
  SUBCASE("test are_points_equal") {
    SUBCASE("test true") {
      math::Vector3 point_1 = {1, 1, 1};
      math::Vector3 point_2 = {1.0, 1.0, 1.0};
      CHECK(math::are_points_equal(point_1, point_2, 0.001) ==
            doctest::Approx(1));
    }
    SUBCASE("test false") {
      math::Vector3 point_1 = {1, 1, 1};
      math::Vector3 point_2 = {-1, -1, 0};
      CHECK(math::are_points_equal(point_1, point_2, 0.001) ==
            doctest::Approx(0));
    }
    SUBCASE("test borderline but outside bounds") {
      math::Vector3 point_1 = {1, 1, 1};
      math::Vector3 point_2 = {1.01, 1.01, 1.01};
      CHECK(math::are_points_equal(point_1, point_2, 0.001) ==
            doctest::Approx(0));
    }
    SUBCASE("test borderline but inside bounds") {
      math::Vector3 point_1 = {1, 1, 1};
      math::Vector3 point_2 = {1.00001, 1.00001, 1.00001};
      CHECK(math::are_points_equal(point_1, point_2, 0.001) ==
            doctest::Approx(1));
    }
  }
  SUBCASE("test are_matrices_equal") {
    SUBCASE("test true") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      CHECK(math::are_matrices_equal(matrix_1, matrix_2) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      math::Matrix3x3 matrix_1 = {1, 2, 1, 2, 1, 2, 2, 2, 2};
      math::Matrix3x3 matrix_2 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
      CHECK(math::are_matrices_equal(matrix_1, matrix_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.01, 1.01, 1.01, 1.01, 1.01,
                                  1.01, 1.01, 1.01, 1.01};
      CHECK(math::are_matrices_equal(matrix_1, matrix_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.00001, 1.00001, 1.00001, 1.00001, 1.00001,
                                  1.00001, 1.00001, 1.00001, 1.00001};
      CHECK(math::are_matrices_equal(matrix_1, matrix_2) == doctest::Approx(1));
    }
  }
  SUBCASE("test roughly_equal") {
    SUBCASE("test true") {
      Real a = 1;
      Real b = 1.0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      Real a = 1;
      Real b = 0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      Real a = 1;
      Real b = 1.01;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      Real a = 1;
      Real b = 1.00001;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
  }
  SUBCASE("test roughly_equal for doubles") {
    SUBCASE("test true") {
      double a = 1;
      double b = 1.0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      double a = 1;
      double b = 0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      double a = 1;
      double b = 1.01;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      double a = 1;
      double b = 1.00001;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
  }
  SUBCASE("test roughly_equal for floats") {
    SUBCASE("test true") {
      float a = 1;
      float b = 1.0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      float a = 1;
      float b = 0;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      float a = 1;
      float b = 1.01;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      float a = 1;
      float b = 1.00001;
      CHECK(math::roughly_equal(a, b, 0.001) == doctest::Approx(1));
    }
  }
  SUBCASE("test roughly_equal for matrices") {
    SUBCASE("test true") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      CHECK(math::roughly_equal(matrix_1, matrix_2) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {2, 2, 2, 2, 2, 2, 2, 2, 2};
      CHECK(math::roughly_equal(matrix_1, matrix_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 1.01};
      CHECK(math::roughly_equal(matrix_1, matrix_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      math::Matrix3x3 matrix_1 = {1, 1, 1, 1, 1, 1, 1, 1, 1};
      math::Matrix3x3 matrix_2 = {1.00001, 1.00001, 1.00001, 1.00001, 1.00001,
                                  1.00001, 1.00001, 1.00001, 1.00001};
      CHECK(math::roughly_equal(matrix_2, matrix_1) == doctest::Approx(1));
    }
  }
  SUBCASE("test roughly_equal for vectors") {
    SUBCASE("test true") {
      math::Vector3 vector_1 = {1, 1, 1};
      math::Vector3 vector_2 = {1.0, 1.0, 1.0};
      CHECK(math::roughly_equal(vector_1, vector_2) == doctest::Approx(1));
    }
    SUBCASE("test false") {
      math::Vector3 vector_1 = {1, 1, 1};
      math::Vector3 vector_2 = {0, 0, 0};
      CHECK(math::roughly_equal(vector_1, vector_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but outside") {
      math::Vector3 vector_1 = {1, 1, 1};
      math::Vector3 vector_2 = {1.01, 1.01, 1.01};
      CHECK(math::roughly_equal(vector_1, vector_2) == doctest::Approx(0));
    }
    SUBCASE("test borderline but inside") {
      math::Vector3 vector_1 = {1, 1, 1};
      math::Vector3 vector_2 = {1.00001, 1.00001, 1.00001};
      CHECK(math::roughly_equal(vector_1, vector_2) == doctest::Approx(1));
    }
  }
}
