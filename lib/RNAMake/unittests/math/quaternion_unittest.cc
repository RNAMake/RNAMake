//
// Created by Joseph Yesselman on 3/4/19.
//


#include "../common.hpp"

#include "base/settings.h"
#include "base/file_io.h"
#include "math/quaternion.h"
#include "math/numerical.h"

typedef std::vector<std::vector<double>> Vector2D;

math::Quaternion
quaternion_from_str(
        String s) {

    auto spl = base::split_str_by_delimiter(s, " ");
    auto values = std::vector<double>();
    for(auto const & e : spl) {
        if(e.size() < 3) { continue; }
        values.push_back(std::stod(e));
    }

    return math::Quaternion(values[0], values[1], values[2], values[3]);
}

TEST_CASE( "Test Quaternion calculations", "[Quaternion]" ) {

    SECTION("test power iteration method for getting eigen values") {
        //auto m = std::vector<std::vector<
        auto m = Vector2D{std::vector<double>{0.5, 0.5},
                          std::vector<double>{0.2, 0.8}};

        SECTION("test dotting with Vector2D") {
            auto v = std::vector<double>{0.5, 0.5};
            auto vr = std::vector<double>(2);
            math::dot_vector(m, v, vr);
            REQUIRE(math::are_floats_equal(vr[0], 0.5));
            REQUIRE(math::are_floats_equal(vr[1], 0.5));
        }

        SECTION("test calculating norm") {
            auto v = std::vector<double>{0.5, 0.5};
            REQUIRE(math::are_floats_equal(math::norm(v), 0.707107));
        }

        SECTION("test full method") {
            auto eigen_values = std::vector<double>(2);
            math::power_iteration(m, eigen_values, 100);
            REQUIRE(math::are_floats_equal(eigen_values[0], 0.707107));
            REQUIRE(math::are_floats_equal(eigen_values[1], 0.707107));
        }
    }

    SECTION("test quaternion averaging") {
        SECTION("test reproducing same quaternion") {
            auto path = base::unittest_resource_dir() + "/math/test_quaternions.dat";
            auto lines = base::get_lines_from_file(path);

            auto averager = math::AverageQuaternionCalculator();
            auto q_init = quaternion_from_str(lines[0]);
            averager.add_quaternion(q_init);

            auto q_avg = averager.get_average();

            for (int i = 0; i < 4; i++) {
                REQUIRE(math::are_floats_equal(q_init[i], q_avg[i]));
            }

            // add 9 more
            for (int i = 0; i < 9; i++) {
                averager.add_quaternion(q_init);
            }

            q_avg = averager.get_average();
            for (int i = 0; i < 4; i++) {
                REQUIRE(math::are_floats_equal(q_init[i], q_avg[i]));
            }
        }

        SECTION("reproduce matlab values") {
            auto path = base::unittest_resource_dir() + "/math/test_quaternions.dat";
            auto lines = base::get_lines_from_file(path);

            auto averager = math::AverageQuaternionCalculator();
            auto q_init = quaternion_from_str(lines[0]);

            for(int i = 1; i < lines.size(); i++) {
                if(lines[i].size() < 10) { continue; }
                auto q = quaternion_from_str(lines[i]);
                averager.add_quaternion(q);
            }

            auto matlab_result_q = math::Quaternion(0.4199, 0.6709, -0.4577, -0.4050);
            auto q_avg = averager.get_average();

            for (int i = 0; i < 4; i++) {
                REQUIRE(math::are_floats_equal(std::abs(matlab_result_q[i]), abs(q_avg[i])));
            }

        }

        SECTION("covert to a rotation matrix and back") {
            //auto q = math::get_random_quaternion();
            //std::cout << q << std::endl;
            auto q = math::Quaternion(0.656091, 0.545157, 0.0359373, 0.520632);
            auto r = q.get_rotation_matrix();

            auto q1 = get_quaternion_from_matrix(r);

            for (int i = 0; i < 4; i++) {
                REQUIRE(math::are_floats_equal(std::abs(q[i]), abs(q1[i])));
            }
        }
    }

}





























