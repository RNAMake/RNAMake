//
// Created by Joe Yesselman on 4/18/22.
//

#include "../common.hpp"

#include <sstream>

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
    }
    SUBCASE("too many input arguments") {
      Reals nums2 = {0.0, 1.0, 2.0, 3.0};
      CHECK_THROWS_AS(math::Vector3{nums2}, base::MathException);
    }
    SUBCASE("too few input arguments"){
      Reals nums3 = {0.0, 1.0};
      Reals nums4 = {1.0};
      CHECK_THROWS_AS(math::Vector3{nums3}, base::MathException);
      CHECK_THROWS_AS(math::Vector3{nums4}, base::MathException);
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
    SUBCASE("using = operator constant") {
      math::Vector3 vec = {0.0f, 1.0f, 2.0f};
      Real new_val = 3;
      vec = new_val;
      CHECK(vec.get_x() == doctest::Approx(3));
      CHECK(vec.get_y() == doctest::Approx(3));
      CHECK(vec.get_z() == doctest::Approx(3));
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
      SUBCASE("const + vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real c = 0.1f;
        vec_1 = c + vec_1;
        CHECK(vec_1.get_x() == doctest::Approx(0.1));
        CHECK(vec_1.get_y() == doctest::Approx(1.1));
        CHECK(vec_1.get_z() == doctest::Approx(2.1));
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
      SUBCASE("constant - vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        vec_1 = 3 - vec_1;
        CHECK(vec_1.get_x() == doctest::Approx(3.0));
        CHECK(vec_1.get_y() == doctest::Approx(2.0));
        CHECK(vec_1.get_z() == doctest::Approx(1.0));
      }
    }
    SUBCASE("test *") {
      SUBCASE("test vec * const") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = vec_1 * 2;
        CHECK(vec_2.get_x() == doctest::Approx(0));
        CHECK(vec_2.get_y() == doctest::Approx(2));
        CHECK(vec_2.get_z() == doctest::Approx(4));
      }
      SUBCASE("test const * vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = 2 * vec_1;
        CHECK(vec_2.get_x() == doctest::Approx(0));
        CHECK(vec_2.get_y() == doctest::Approx(2));
        CHECK(vec_2.get_z() == doctest::Approx(4));
      }
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
    SUBCASE("test -=") {
      SUBCASE("test vec -= vec") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        vec_1 -= vec_2;
        CHECK(vec_1.get_x() == doctest::Approx(0));
        CHECK(vec_1.get_y() == doctest::Approx(0));
        CHECK(vec_1.get_z() == doctest::Approx(0));
      }
      SUBCASE("test vec -= const") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        vec_1 -= 1.0f;
        CHECK(vec_1.get_x() == doctest::Approx(-1));
        CHECK(vec_1.get_y() == doctest::Approx(0));
        CHECK(vec_1.get_z() == doctest::Approx(1));
      }
    }
    SUBCASE("test *= const") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      vec_1 *= 2;
      CHECK(vec_1.get_x() == doctest::Approx(0));
      CHECK(vec_1.get_y() == doctest::Approx(2));
      CHECK(vec_1.get_z() == doctest::Approx(4));
    }
    SUBCASE("test /= const") {
      SUBCASE("trival test") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        vec_1 /= 2;
        CHECK(vec_1.get_x() == doctest::Approx(0));
        CHECK(vec_1.get_y() == doctest::Approx(0.5f));
        CHECK(vec_1.get_z() == doctest::Approx(1.0f));
      }
      SUBCASE("divide by zero test") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        CHECK_THROWS_AS(vec_1 /= 0, base::MathException);
      }
    }
    SUBCASE("test == ") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
      CHECK(vec_1 == vec_2);
    }
    SUBCASE("test != ") {
      SUBCASE("trival test") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.1f, 1.0f, 2.0f};
        CHECK(vec_1 != vec_2);
      }
      SUBCASE("really small change") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.000001f, 2.0f};
        CHECK(vec_1 != vec_2);
      }
    }
    SUBCASE("test <<") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      std::stringstream ss;
      ss << vec_1;
      CHECK(ss.str() == "0.000000 1.000000 2.000000");
    }
    SUBCASE("test cross") {
      SUBCASE("trivial test") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        math::Vector3 cross = vec_1.cross(vec_2);
        CHECK(cross.get_x() == doctest::Approx(0));
        CHECK(cross.get_y() == doctest::Approx(0));
        CHECK(cross.get_z() == doctest::Approx(0));
      }
      SUBCASE("more complex test") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {1.0f, 1.0f, 1.0f};
        math::Vector3 cross = vec_1.cross(vec_2);
        CHECK(cross.get_x() == doctest::Approx(-1));
        CHECK(cross.get_y() == doctest::Approx(2));
        CHECK(cross.get_z() == doctest::Approx(-1));
      }
      SUBCASE("slightly more complex test") {
        math::Vector3 vec_1 = {3.0f, -5.0f, 1.0f};
        math::Vector3 vec_2 = {6.0f, 4.0f, -1.0f};
        math::Vector3 cross = vec_1.cross(vec_2);
        CHECK(cross.get_x() == doctest::Approx(1));
        CHECK(cross.get_y() == doctest::Approx(9));
        CHECK(cross.get_z() == doctest::Approx(42));
      }
    }
    SUBCASE("test negating") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_1_negate = vec_1.negate();
      CHECK(vec_1_negate.get_x() == doctest::Approx(0));
      CHECK(vec_1_negate.get_y() == doctest::Approx(-1));
      CHECK(vec_1_negate.get_z() == doctest::Approx(-2));
    }
    SUBCASE("test negating and copying") {
      math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
      math::Vector3 vec_2 = vec_1.negated();
      CHECK(vec_2.get_x() == doctest::Approx(0));
      CHECK(vec_2.get_y() == doctest::Approx(-1));
      CHECK(vec_2.get_z() == doctest::Approx(-2));
    }
    SUBCASE("test distance function") {
      math::Vector3 vec_1 = {0.0f, 0.0f, 0.0f};
      math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
      Real vec_2_distance = vec_1.distance(vec_2);
      CHECK(vec_2_distance == doctest::Approx(sqrt(5)));
    }
    SUBCASE("test distance sqaured function") {
      math::Vector3 vec_1 = {0.0f, 0.0f, 0.0f};
      math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
      Real vec_2_distance_squared = vec_1.distance_squared(vec_2);
      CHECK(vec_2_distance_squared == doctest::Approx(5));
    }
    SUBCASE("test dot products"){
      SUBCASE("test zero vector dot products") {
        math::Vector3 vec_1 = {};
        math::Vector3 vec_2 = {};
        Real dot_product = vec_1.dot(vec_2);
        CHECK(dot_product == doctest::Approx(0));
      }
      SUBCASE("test simple dot product") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        math::Vector3 vec_2 = {0.0f, 1.0f, 2.0f};
        Real dot_product = vec_1.dot(vec_2);
        CHECK(dot_product == doctest::Approx(5));
      }
      SUBCASE("test complex dot product") {
        math::Vector3 vec_1 = {-1.0f, 2.0f, 4.0f};
        math::Vector3 vec_2 = {6.0f, -3.0f, 9.0f};
        Real dot_product = vec_1.dot(vec_2);
        CHECK(dot_product == doctest::Approx(24));
      }
    }
    SUBCASE("test normalization"){
      SUBCASE("test zero vector") {
        math::Vector3 vec_1 = {};
        math::Vector3 vec_2 = {0.0f, 0.0f, 0.0f};
        CHECK_THROWS_AS(vec_1.normalize(), base::MathException);
        CHECK_THROWS_AS(vec_2.normalize(), base::MathException);
      }
      SUBCASE("test simple vector") {
        math::Vector3 vec_1 = {1.0f, 2.0f, 3.0f};
        math::Vector3 vec_1_normalized = vec_1.normalize();
        CHECK(vec_1_normalized.get_x() == doctest::Approx(1/sqrt(14)));
        CHECK(vec_1_normalized.get_y() == doctest::Approx(2/sqrt(14)));
        CHECK(vec_1_normalized.get_z() == doctest::Approx(3/sqrt(14)));
      }
      SUBCASE("test complex vector") {
        math::Vector3 vec_1 = {3.0f, -5.0f, 4.0f};
        math::Vector3 vec_1_normalized = vec_1.normalize();
        CHECK(vec_1_normalized.get_x() == doctest::Approx(3/sqrt(50)));
        CHECK(vec_1_normalized.get_y() == doctest::Approx(-5/sqrt(50)));
        CHECK(vec_1_normalized.get_z() == doctest::Approx(4/sqrt(50)));
      }
    }
    SUBCASE("test zeroing functions") {
      SUBCASE("test zero vectors") {
        math::Vector3 vec_1 = {};
        math::Vector3 vec_2 = {0.0f, 0.0f, 0.0f};
        math::Vector3 vec_1_zeroed = vec_1.zero();
        math::Vector3 vec_2_zeroed = vec_2.zero();
        CHECK(vec_1_zeroed.get_x() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_y() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_z() == doctest::Approx(0));
        CHECK(vec_2_zeroed.get_x() == doctest::Approx(0));
        CHECK(vec_2_zeroed.get_y() == doctest::Approx(0));
        CHECK(vec_2_zeroed.get_z() == doctest::Approx(0));
      }
      SUBCASE("test tirvial vector") {
        math::Vector3 vec_1 = {1.0f, 2.0f, 3.0f};
        math::Vector3 vec_1_zeroed = vec_1.zero();
        CHECK(vec_1_zeroed.get_x() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_y() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_z() == doctest::Approx(0));
      }
      SUBCASE("test complex vector") {
        math::Vector3 vec_1 = {3.0f, -5.0f, 2.0f};
        math::Vector3 vec_1_zeroed = vec_1.zero();
        CHECK(vec_1_zeroed.get_x() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_y() == doctest::Approx(0));
        CHECK(vec_1_zeroed.get_z() == doctest::Approx(0));
      }
    }
    SUBCASE("test getting functions") {
      SUBCASE("test getting length") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_length = vec_1.get_length();
        CHECK(vec_1_length == doctest::Approx(sqrt(5)));
      }
      SUBCASE("test getting length squared") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_length_squared = vec_1.get_length_squared();
        CHECK(vec_1_length_squared == doctest::Approx(5));
      }
      SUBCASE("test get norm") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_norm = vec_1.get_norm();
        CHECK(vec_1_norm == doctest::Approx(sqrt(5)));
      }
      SUBCASE("test getting norm squared") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_norm_squared = vec_1.get_norm_squared();
        CHECK(vec_1_norm_squared == doctest::Approx(5));
      }
      SUBCASE("test get magnitude") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_magnitude = vec_1.get_magnitude();
        CHECK(vec_1_magnitude == doctest::Approx(sqrt(5)));
      }
      SUBCASE("test getting magnitude squared") {
        math::Vector3 vec_1 = {0.0f, 1.0f, 2.0f};
        Real vec_1_magnitude_squared = vec_1.get_magnitude_squared();
        CHECK(vec_1_magnitude_squared == doctest::Approx(5));
      }
      SUBCASE("test getting vector components") {
        math::Vector3 vec_1 = {1.0f, 2.0f, 3.0f};
        Real vec_1_x = vec_1.get_x();
        Real vec_1_y = vec_1.get_y();
        Real vec_1_z = vec_1.get_z();
        CHECK(vec_1_x == doctest::Approx(1));
        CHECK(vec_1_y == doctest::Approx(2));
        CHECK(vec_1_z == doctest::Approx(3));
      }
    }
    SUBCASE("test setting fxns") {
      SUBCASE("set vector components to zero") {
        math::Vector3 vec_1 = {134.4f, 451.2f, 23.2f};
        vec_1.set_x(0);
        vec_1.set_y(0);
        vec_1.set_z(0);
        CHECK(vec_1.get_x() == doctest::Approx(0));
        CHECK(vec_1.get_y() == doctest::Approx(0));
        CHECK(vec_1.get_z() == doctest::Approx(0));
      }
      SUBCASE("assign values to empty vector") {
        math::Vector3 vec_1 = {};
        vec_1.set_x(1);
        vec_1.set_y(2);
        vec_1.set_z(3);
        CHECK(vec_1.get_x() == doctest::Approx(1));
        CHECK(vec_1.get_y() == doctest::Approx(2));
        CHECK(vec_1.get_z() == doctest::Approx(3));
      }
    }
    SUBCASE("test vector from string") {
      SUBCASE("test normal circumstance") {
        String vec_1_string = "1.0 2.0 3.0";
        math::Vector3 vec_1 = math::vector_from_str(vec_1_string);
        CHECK(vec_1.get_x() == doctest::Approx(1));
        CHECK(vec_1.get_y() == doctest::Approx(2));
        CHECK(vec_1.get_z() == doctest::Approx(3));
      }
      SUBCASE("test too few inputs") {
        String vec_1_string = "1.0 2.0";
        CHECK_THROWS_AS(math::vector_from_str(vec_1_string), base::InputException);
      }
      SUBCASE("test too many inputs") {
        String vec_1_string = "1.0 2.0 3.0 4.0";
        CHECK_THROWS_AS(math::vector_from_str(vec_1_string), base::InputException);
      }
      SUBCASE("tset zero input") {
        String vec_1_string = "";
        CHECK_THROWS_AS(math::vector_from_str(vec_1_string), base::InputException);
      }
    }
    SUBCASE("test constructing vectors string") {
      SUBCASE("test normal circumstance") {
        String vec_1_string = "1 2 3";
        math::Vector3s vecs;
        math::vectors_from_str(vec_1_string, vecs);
        CHECK(vecs.size() == 1);
        CHECK(vecs[0].get_x() == doctest::Approx(1));
        CHECK(vecs[0].get_y() == doctest::Approx(2));
        CHECK(vecs[0].get_z() == doctest::Approx(3));
      }
      SUBCASE("test not enough inputs") {
        String vec_1_string = "1 2";
        math::Vector3s vec_1 = {};
        CHECK_THROWS_AS(math::vectors_from_str(vec_1_string, vec_1), base::InputException);
      }
      SUBCASE("test too many inputs") {
        String vec_1_string = "1 2 3 4";
        math::Vector3s vec_1 = {};
        CHECK_THROWS_AS(math::vectors_from_str(vec_1_string, vec_1), base::InputException);
      }
      SUBCASE("test no inputs") {
        String vec_1_string;
        math::Vector3s vec_1 = {};
        CHECK_THROWS_AS(math::vectors_from_str(vec_1_string, vec_1), base::InputException);
      }
    }
  }
}
