

#include "../common.hpp"

#include "structure/atom.h"
#include "structure/is_equal.h"


TEST_CASE( "Test Atoms for Structure") {
    auto a = std::make_shared<structure::Atom>("P", math::Point(0, 1, 2));
    auto p = math::Point(0, 1, 2);
    
    CHECK(a->coords() == p);
    CHECK(a->name() == "P");
    
    SUBCASE("do atoms have correct pdb output") {
        auto s = a->to_pdb_str(1);
        auto ref = "ATOM      1  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
        
        CHECK(s == ref);
        
        auto s1 = a->to_pdb_str(10);
        auto ref1 = "ATOM     10  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
        
        CHECK(s1 == ref1);
    }
    
    
    SUBCASE("do atoms stringify properly") {
        auto s = a->to_str();
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

