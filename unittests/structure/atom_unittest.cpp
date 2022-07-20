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

    CHECK((atom_1.get_name() == atom_6.get_name()) == true);
    CHECK((atom_1.get_name() == atom_3.get_name()) == false);

    CHECK(atom_1.get_x() == doctest::Approx(0));
    CHECK(atom_1.get_y() == doctest::Approx(0));
    CHECK(atom_1.get_z() == doctest::Approx(0));

    CHECK(atom_4.get_x() == doctest::Approx(-10));
    CHECK(atom_4.get_y() == doctest::Approx(10));
    CHECK(atom_4.get_z() == doctest::Approx(10));

    CHECK((atom_1.get_name() == atom_6.get_name()) == true);
    CHECK((atom_2.get_name() == atom_3.get_name()) == false);
  }
  SUBCASE("test atom construction error") {
    math::Vector3 position = math::Vector3(0, 0, 0);
    String name = "";
    //CHECK_THROWS_AS(Atom(name, position), base::InputException);
  }
}