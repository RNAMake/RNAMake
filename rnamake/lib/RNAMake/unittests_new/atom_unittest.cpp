#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "is_equal.hpp"

#include "structure/atom.h"


TEST_CASE( "Test Atoms for Structure", "[Atom]" ) {
    auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
    auto p = Point(0, 1, 2);
    
    REQUIRE(a->coords() == p);
    REQUIRE(a->name() == "P");
    
    SECTION("do atoms have correct pdb output") {
        auto s = a->to_pdb_str(1);
        auto ref = "ATOM      1  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
        
        REQUIRE(s == ref);
        
        auto s1 = a->to_pdb_str(10);
        auto ref1 = "ATOM     10  P   C   A   1       0.000   1.000   2.000  1.00 62.18           P\n";
        
        REQUIRE(s1 == ref1);
    }
    
    
    SECTION("do atoms stringify properly") {
        auto s = a->to_str();
        auto a2 = std::make_shared<Atom>(s);
        
        REQUIRE(s == "P 0 1 2");
        REQUIRE(are_atoms_equal(a, a2));
        
    }
    
    SECTION("do atoms copy properly") {
        auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
        auto a2 = std::make_shared<Atom>(*a);

        REQUIRE(are_atoms_equal(a, a2));
    }
    
}




