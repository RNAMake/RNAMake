

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/residue.h"

TEST_CASE( "Test Residues for Secondary Structure" ) {
    
    auto ss_r = std::make_shared<secondary_structure::Residue>("G", ".", 1, "A", util::Uuid());
    
    CHECK(ss_r->res_type() == secondary_structure::ResType::GUA);
    
    SUBCASE("check for valid residue name") {
        REQUIRE_THROWS_AS(secondary_structure::Residue("K", ".", 1, "A", util::Uuid()),
                          secondary_structure::Exception);
    }
    
    SUBCASE("test residue copy") {
        auto ss_r_copy = std::make_shared<secondary_structure::Residue>(*ss_r);
        
        CHECK(ss_r_copy->name() == ss_r->name());
        CHECK(ss_r_copy->res_type() == ss_r->res_type());

    }
    
    SUBCASE("test stringifying residue" ) {
        auto s = ss_r->to_str();
        auto ss_r2 = std::make_shared<secondary_structure::Residue>(s);
        
        CHECK(ss_r2->name() == ss_r->name());
        CHECK(ss_r2->res_type() == ss_r->res_type());

        SUBCASE("catch bad string to build from") {
            REQUIRE_THROWS_AS(secondary_structure::Residue(""), secondary_structure::Exception);
        }
    }

}