

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/structure.h"

TEST_CASE( "Test Structure for Secondary Structure" ) {
    
    SUBCASE("test creation of structure from sequence and structure") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        
        CHECK(s.sequence() == "GGACC");
        CHECK(s.dot_bracket() == "((.))");
        
        SUBCASE("sequence and structure not the same length") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACCC", "((.))"),
                              secondary_structure::Exception);
        }
        
        SUBCASE("invalid structure start character") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACC", ")(.))"),
                              secondary_structure::Exception);

        }
        
        SUBCASE("invalid residue name in sequence") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGKCC", "((.))"),
                              secondary_structure::Exception);
            
        }
        
        SUBCASE("invalid residue structure element in structure") {
            REQUIRE_THROWS_AS(secondary_structure::Structure("GGACC", "((*))"),
                              secondary_structure::Exception);
            
        }
        

        SUBCASE("has to recycle chain ids since there are more than 26 chains") {
            auto seq = String(), ss = String();
            for(int i = 0; i < 100; i++) {
                seq += "A&";
                ss += ".&";
            }
        
            auto s1 = secondary_structure::Structure(seq, ss);
            CHECK(s1.chains().size() == 100);
            CHECK(s1.chains()[99]->residues()[0]->chain_id() == "D");
        }
        
    }
    
    SUBCASE("test copy constructor") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto s_copy = secondary_structure::Structure(s);
        
        CHECK(s_copy.sequence() == "GGACC");
        CHECK(s_copy.dot_bracket() == "((.))");

    }
    
    SUBCASE("test stringify of structure") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto str = s.to_str();
        auto s_copy = secondary_structure::Structure(str);
        
        CHECK(s_copy.sequence() == "GGACC");
        CHECK(s_copy.dot_bracket() == "((.))");

    }
    
    SUBCASE("test ability to find residue") {
        auto s = secondary_structure::Structure("GGACC", "((.))");
        auto r = s.get_residue(1, "A", "");
        CHECK(r != nullptr);
        
        auto r2 = s.get_residue(r->uuid());
        CHECK(r != nullptr);
        
        CHECK(s.get_residue(99, "A", "") == nullptr);
        CHECK(s.get_residue(1, "B", "") == nullptr);
        CHECK(s.get_residue(util::Uuid()) == nullptr);
        
        
        
    }
    
}