//
// Created by Joe Yesselman on 7/20/22.
//

#include "../common.hpp"
#include <structure/all_atom/atom.h>


TEST_CASE("test atom functions ") {

  using namespace structure::all_atom;

  SUBCASE("test atom construction") {
    math::Vector3 position = math::Vector3(0, 0, 0);
    String name = "hydrogen";
    Atom atom_1 = Atom(name, position);

    math::Vector3 position_2 = math::Vector3(0, 0, 0);
    String name_2 = "helium";
    Atom atom_2 = Atom(name_2, position_2);

    math::Vector3 position_3 = math::Vector3(5, -5, 5);
    String name_3 = "lithium";
    Atom atom_3 = Atom(name_3, position_3);

    math::Vector3 position_4 = math::Vector3(-10, 10, 10);
    String name_4 = "beryllium";
    Atom atom_4 = Atom(name_4, position_4);

    math::Vector3 position_5 = math::Vector3(25, 25, -25);
    String name_5 = "boron";
    Atom atom_5 = Atom(name_5, position_5);

    math::Vector3 position_6 = math::Vector3(0, 0, 0);
    String name_6 = "hydrogen";
    Atom atom_6 = Atom(name_6, position_6);

    math::Vector3 position_7 = math::Vector3(0.001, 0.001, 0.001);
    String name_7 = "carbon";
    Atom atom_7 = Atom(name_7, position_7);

    math::Vector3 position_8 = math::Vector3(0, 0, 0);
    String name_8 = "HYDROGEN";
    Atom atom_8 = Atom(name_8, position_8);

    // checks equality fxns and some basic logic
    CHECK((atom_1.get_name() == atom_6.get_name()) == true);
    CHECK((atom_1.get_name() == atom_3.get_name()) == false);
    // TODO should names be transferred to lowercase?
    CHECK((atom_1.get_name() == atom_8.get_name()) == false);

    CHECK((atom_1.get_coords() == atom_6.get_coords()) == true);
    CHECK((atom_1.get_coords() == atom_7.get_coords()) == false);
    CHECK((atom_1.get_coords() == atom_3.get_coords()) == false);

    // checks atom coord-getting fxn
    CHECK(atom_1.get_x() == doctest::Approx(0));
    CHECK(atom_1.get_y() == doctest::Approx(0));
    CHECK(atom_1.get_z() == doctest::Approx(0));

    CHECK(atom_2.get_x() == doctest::Approx(0));
    CHECK(atom_2.get_y() == doctest::Approx(0));
    CHECK(atom_2.get_z() == doctest::Approx(0));

    CHECK(atom_3.get_x() == doctest::Approx(5));
    CHECK(atom_3.get_y() == doctest::Approx(-5));
    CHECK(atom_3.get_z() == doctest::Approx(5));

    CHECK(atom_4.get_x() == doctest::Approx(-10));
    CHECK(atom_4.get_y() == doctest::Approx(10));
    CHECK(atom_4.get_z() == doctest::Approx(10));

    CHECK(atom_5.get_x() == doctest::Approx(25));
    CHECK(atom_5.get_y() == doctest::Approx(25));
    CHECK(atom_5.get_z() == doctest::Approx(-25));

    CHECK(atom_6.get_x() == doctest::Approx(0));
    CHECK(atom_6.get_y() == doctest::Approx(0));
    CHECK(atom_6.get_z() == doctest::Approx(0));

    CHECK(atom_7.get_x() == doctest::Approx(0.001));
    CHECK(atom_7.get_y() == doctest::Approx(0.001));
    CHECK(atom_7.get_z() == doctest::Approx(0.001));

    // checks get_coords functions
    CHECK(atom_1.get_coords() == position);
    CHECK(atom_2.get_coords() == position_2);
    CHECK(atom_3.get_coords() == position_3);
    CHECK(atom_4.get_coords() == position_4);
    CHECK(atom_5.get_coords() == position_5);
    CHECK(atom_6.get_coords() == position_6);
    CHECK(atom_7.get_coords() == position_7);

    // double-checks get_coords functions
    math::Vector3 vector_1 = atom_1.get_coords();
    math::Vector3 vector_2 = atom_2.get_coords();
    math::Vector3 vector_3 = atom_3.get_coords();
    math::Vector3 vector_4 = atom_4.get_coords();
    math::Vector3 vector_5 = atom_5.get_coords();
    math::Vector3 vector_6 = atom_6.get_coords();
    math::Vector3 vector_7 = atom_7.get_coords();

    CHECK(vector_1.get_x() == doctest::Approx(0));
    CHECK(vector_1.get_y() == doctest::Approx(0));
    CHECK(vector_1.get_z() == doctest::Approx(0));

    CHECK(vector_2.get_x() == doctest::Approx(0));
    CHECK(vector_2.get_y() == doctest::Approx(0));
    CHECK(vector_2.get_z() == doctest::Approx(0));

    CHECK(vector_3.get_x() == doctest::Approx(5));
    CHECK(vector_3.get_y() == doctest::Approx(-5));
    CHECK(vector_3.get_z() == doctest::Approx(5));

    CHECK(vector_4.get_x() == doctest::Approx(-10));
    CHECK(vector_4.get_y() == doctest::Approx(10));
    CHECK(vector_4.get_z() == doctest::Approx(10));

    CHECK(vector_5.get_x() == doctest::Approx(25));
    CHECK(vector_5.get_y() == doctest::Approx(25));
    CHECK(vector_5.get_z() == doctest::Approx(-25));

    CHECK(vector_6.get_x() == doctest::Approx(0));
    CHECK(vector_6.get_y() == doctest::Approx(0));
    CHECK(vector_6.get_z() == doctest::Approx(0));

    CHECK(vector_7.get_x() == doctest::Approx(0.001));
    CHECK(vector_7.get_y() == doctest::Approx(0.001));
    CHECK(vector_7.get_z() == doctest::Approx(0.001));

  }
  /*
  SUBCASE("test atom construction error") {
    math::Vector3 position = math::Vector3(0, 0, 0);
    String name = "";
    // CHECK_THROWS_AS(Atom(name, position), base::InputException);
  }
   */
}