//
// Created by Joseph Yesselman on 2019-04-18.
//

#include "../common.hpp"

#include <math/xyz_matrix.h>
#include <math/euler.h>

math::Matrix
rotation_about_x_axis(
        float degrees) {
    float a = M_PI/180*degrees;
    return math::Matrix(1.0, 0.0, 0.0,
                        0.0, cos(a), sin(a),
                        0.0, -sin(a), cos(a));
}


TEST_CASE( "Test Rotations ", "[Rotations]" ) {

    SECTION("test axis angle rotations") {
        auto m = math::Matrix( 0.866, 0.5, 0,
                              -0.5, 0.866, 0,
                              0, 0, 1);
        auto aa = math::AxisAngle();
        math::axis_angle_from_matrix(m, aa);

        REQUIRE(abs(30 - math::degrees(aa.angle)) < 0.01);
        REQUIRE(aa.axis.distance(math::Point(0, 0, -1)) < 0.01);

        m = math::Matrix( 0.866, -0.5, 0,
                          0.5, 0.866, 0,
                          0, 0, 1);

        math::axis_angle_from_matrix(m, aa);
        REQUIRE(abs(30 - math::degrees(aa.angle)) < 0.01);
        REQUIRE(aa.axis.distance(math::Point(0, 0, 1)) < 0.01);

    }


}