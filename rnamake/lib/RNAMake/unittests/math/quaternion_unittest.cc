//
// Created by Joseph Yesselman on 3/4/19.
//


#include "../common.hpp"

#include "base/settings.h"
#include "base/file_io.h"
#include "math/quaternion.h"
#include "math/numerical.h"

typedef std::vector<std::vector<double>> Vector2D;

TEST_CASE( "Test Quaternion calculations", "[Quaternion]" ) {

    SECTION("test power iteration method for getting eigen values") {
        //auto m = std::vector<std::vector<
        auto m = Vector2D{std::vector<double>{0.5, 0.5},
                          std::vector<double>{0.2, 0.8}};

        SECTION("test dotting with Vector2D") {
            auto v = std::vector<double>{0.5, 0.5};
            auto vr = std::vector<double>(2);
            dot_vector(m, v, vr);
            REQUIRE(are_floats_equal(vr[0], 0.5));
            REQUIRE(are_floats_equal(vr[1], 0.5));
        }

        SECTION("test calculating norm") {
            auto v = std::vector<double>{0.5, 0.5};
            REQUIRE(are_floats_equal(norm(v), 0.707107));
        }

        SECTION("test full method") {
            auto eigen_values = std::vector<double>(2);
            power_iteration(m, eigen_values, 100);
            REQUIRE(are_floats_equal(eigen_values[0], 0.707107));
            REQUIRE(are_floats_equal(eigen_values[1], 0.707107));
        }
    }


}