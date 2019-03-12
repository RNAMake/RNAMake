


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/residue.h"

TEST_CASE( "Test Residues for Secondary Structure", "[SSResidue]" ) {
    
    auto ss_r = std::make_shared<secondary_structure::Residue>("G", ".", 1, "A", util::Uuid());
    
    REQUIRE(ss_r->res_type() == 2);
    
    SECTION("check for valid residue name") {
        REQUIRE_THROWS_AS(secondary_structure::Residue("K", ".", 1, "A", util::Uuid()),
                          secondary_structure::Exception);
    }
    
    SECTION("test residue copy") {
        auto ss_r_copy = std::make_shared<secondary_structure::Residue>(*ss_r);
        
        REQUIRE(ss_r_copy->name() == ss_r->name());
        REQUIRE(ss_r_copy->res_type() == ss_r->res_type());

    }
    
    SECTION("test stringifying residue" ) {
        auto s = ss_r->to_str();
        auto ss_r2 = std::make_shared<secondary_structure::Residue>(s);
        
        REQUIRE(ss_r2->name() == ss_r->name());
        REQUIRE(ss_r2->res_type() == ss_r->res_type());

        SECTION("catch bad string to build from") {
            REQUIRE_THROWS_AS(secondary_structure::Residue(""), secondary_structure::Exception);
        }
    }

}