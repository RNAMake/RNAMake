#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include "structure/atom.h"


TEST_CASE( "Test Atoms for Structure", "[Atom]" ) {
    auto a = std::make_shared<Atom>("P", Point(0, 1, 2));
    auto p = Point(0, 1, 2);
    
    REQUIRE(a->coords() == p);
    REQUIRE(a->name() == "P");
    
    SECTION("producing the correct pdb output") {
        auto s = a->to_pdb_str(1);
        auto ref = "ATOM      1  P   C   A   1       0.000   1.000   3.000  1.00 62.18           P\n";
        
        REQUIRE(s == ref);
    }
    
    
}




