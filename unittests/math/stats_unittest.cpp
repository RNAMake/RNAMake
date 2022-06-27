//
// Created by Joe Yesselman on 6/24/22.
//

#include "../common.hpp"
#include <math/stats.hpp>
#include <math/vector_3.hpp>

TEST_CASE("test stat functions ") {
  SUBCASE("test sum function") {
    SUBCASE("test simple") {
      std::vector<double> vector_1 = {1, 2, 3};
      double sum = math::sum(vector_1);
      CHECK(sum == doctest::Approx(6));
    }
    SUBCASE("test complex") {
      std::vector<double> vector_1 = {0, -4.2, 6.0};
      double sum = math::sum(vector_1);
      CHECK(sum == doctest::Approx(1.8));
    }
  }
  SUBCASE("test sqsum function") {
    SUBCASE("test simple") {
      std::vector<double> vector_1 = {1, 2, 3};
      double sum = math::sqsum(vector_1);
      CHECK(sum == doctest::Approx(14));
    }
    SUBCASE("test complex") {
      std::vector<double> vector_1 = {0, -4.2, 6.0};
      double sum = math::sqsum(vector_1);
      CHECK(sum == doctest::Approx(53.64));
    }
  }
  SUBCASE("test standard deviation") {
    SUBCASE("test simple") {
      std::vector<double> vector_1 = {0, 1, 2, 3, 4, 5};
      double sum = math::stdev(vector_1);
      CHECK(sum == doctest::Approx(1.7078));
    }
    SUBCASE("test complex") {
      std::vector<double> vector_1 = {1.7, -9.4, 8.3, 17.9, 2.1, 0.555};
      double sum = math::stdev(vector_1);
      CHECK(sum == doctest::Approx(8.27324));
    }
  }
  SUBCASE("test mean function") {
    SUBCASE("test simple") {
      std::vector<double> vector_1 = {0, 1, 2, 3, 4, 5};
      double sum = math::mean(vector_1);
      CHECK(sum == doctest::Approx(2.5));
    }
    SUBCASE("test complex") {
      std::vector<double> vector_1 = {1.7, -9.4, 8.3, 17.9, 2.1, 0.555};
      double sum = math::mean(vector_1);
      CHECK(sum == doctest::Approx(3.52583));
    }
  }
  SUBCASE("pearson coefficient function") {
    SUBCASE("test simple") {
      std::vector<double> x_axis = {0, 1, 2, 3, 4, 5};
      std::vector<double> y_axis = {0, 1, 2, 3, 4, 5};
      double pearson = math::pearson_coeff(x_axis, y_axis);
      CHECK(pearson == doctest::Approx(1));
    }
    SUBCASE("test complex") {
      std::vector<double> x_axis = {1.7, -9.4, 8.3, 17.9, 2.1, 0.555};
      std::vector<double> y_axis = {2.3, 8.44, 0, -2, 9.21, 3.1};
      double pearson = math::pearson_coeff(x_axis, y_axis);
      CHECK(pearson == doctest::Approx(-0.80671));
    }
  }
  SUBCASE("test avg abval diff") {
    SUBCASE("test super simple") {
      std::vector<double> x_axis = {1, 2, 3, 4, 5, 6};
      std::vector<double> y_axis = {1, 2, 3, 4, 5, 6};
      double sum = math::avg_unsigned_diff(x_axis, y_axis);
      CHECK(sum == doctest::Approx(0));
    }
    SUBCASE("test simple") {
      std::vector<double> x_axis = {1, 2, 3, 4, 5, 6};
      std::vector<double> y_axis = {-1, -2, -3, -4, -5, -6};
      double sum = math::avg_unsigned_diff(x_axis, y_axis);
      CHECK(sum == doctest::Approx(7));
    }
    SUBCASE("test complex") {
      std::vector<double> x_axis = {3.6, 4.7, 8.9, 1.4, 2.3};
      std::vector<double> y_axis = {4.7, 5.6, 6.0, 7.9, 9.1};
      double sum = math::avg_unsigned_diff(x_axis, y_axis);
      CHECK(sum == doctest::Approx(3.64));
    }
  }

}