


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/structure.h"

TEST_CASE( "Test Structure for Secondary Structure", "[SSStructure]" ) {
    
    SECTION("test creation of structure from sequence and structure") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        
        REQUIRE(s.sequence() == "GGACC");
        REQUIRE(s.dot_bracket() == "((.))");
        
        SECTION("sequence and structure not the same length") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACCC", "((.))"),
                              secondary_structure::Exception);
        }
        
        SECTION("invalid structure start character") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACC", ")(.))"),
                              secondary_structure::Exception);

        }
        
        SECTION("invalid residue name in sequence") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGKCC", "((.))"),
                              secondary_structure::Exception);
            
        }
        
        SECTION("invalid residue structure element in structure") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACC", "((*))"),
                              secondary_structure::Exception);
            
        }
        

        SECTION("has to recycle chain ids since there are more than 26 chains") {
            auto seq = String(), ss = String();
            for(int i = 0; i < 100; i++) {
                seq += "A&";
                ss += ".&";
            }
        
            auto s1 = secondary_structure::Structure(seq, ss);
            REQUIRE(s1.chains().size() == 100);
            REQUIRE(s1.chains()[99]->residues()[0]->chain_id() == "D");
        }
        
    }
    
    SECTION("test copy constructor") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto s_copy = secondary_structure::Structure(s);
        
        REQUIRE(s_copy.sequence() == "GGACC");
        REQUIRE(s_copy.dot_bracket() == "((.))");

    }
    
    SECTION("test stringify of structure") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto str = s.to_str();
        auto s_copy = secondary_structure::Structure(str);
        
        REQUIRE(s_copy.sequence() == "GGACC");
        REQUIRE(s_copy.dot_bracket() == "((.))");

    }
    
    SECTION("test ability to find residue") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto r = s.get_residue(1, "A", "");
        REQUIRE(r != nullptr);
        
        auto r2 = s.get_residue(r->uuid());
        REQUIRE(r != nullptr);
        
        REQUIRE(s.get_residue(99, "A", "") == nullptr);
        REQUIRE(s.get_residue(1, "B", "") == nullptr);
        REQUIRE(s.get_residue(util::Uuid()) == nullptr);
        
        
        
    }
    
}