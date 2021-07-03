

#include "../common.hpp"

#include "structure/atom.h"
#include "structure/is_equal.h"


TEST_CASE( "Test Atoms for Structure") {
//TODO Uncomment these after fixing xyz_vector and xyz_matrix

    auto a = std::make_shared<structure::Atom>("P", math::Point(0, 1, 2));
    auto p = math::Point(0, 1, 2);

    CHECK(a->get_coords() == p);
    CHECK(a->get_name() == "P");
//
//    SUBCASE("do atoms have correct pdb output") {
//        auto s = a->get_pdb_str(1);
//        auto ref = "ATOM      1  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
//
//        CHECK(s == ref);
//
//        auto s1 = a->get_pdb_str(10);
//        auto ref1 = "ATOM     10  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
//
//        CHECK(s1 == ref1);
//    }
    
    
    SUBCASE("do atoms stringify properly") {
        auto s = a->get_str();
        auto a2 = std::make_shared<structure::Atom>(s);

        CHECK(s == "P 0 1 2");
        CHECK(are_atoms_equal(a, a2));

    }
    
    SUBCASE("do atoms copy properly") {
        auto a = std::make_shared<structure::Atom>("P", math::Point(0, 1, 2));
        auto a2 = std::make_shared<structure::Atom>(*a);

        CHECK(are_atoms_equal(a, a2));
    }
    
}

