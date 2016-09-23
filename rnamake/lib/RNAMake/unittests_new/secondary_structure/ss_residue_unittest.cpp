


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/residue.h"

TEST_CASE( "Test Residues for Secondary Structure", "[SSResidue]" ) {
    
    auto ss_r = std::make_shared<sstruct::Residue>("G", ".", 1, "A", Uuid());
    
    REQUIRE(ss_r->res_type() == 2);
    
    SECTION("check for valid residue name") {
        REQUIRE_THROWS_AS(sstruct::Residue("K", ".", 1, "A", Uuid()),
                          sstruct::SecondaryStructureException);
    }
    
    SECTION("test residue copy") {
        auto ss_r_copy = std::make_shared<sstruct::Residue>(*ss_r);
        
        REQUIRE(ss_r_copy->name() == ss_r->name());
        REQUIRE(ss_r_copy->res_type() == ss_r->res_type());

    }
    
    SECTION("test stringifying residue" ) {
        auto s = ss_r->to_str();
        auto ss_r2 = std::make_shared<sstruct::Residue>(s);
        
        REQUIRE(ss_r2->name() == ss_r->name());
        REQUIRE(ss_r2->res_type() == ss_r->res_type());

        SECTION("catch bad string to build from") {
            REQUIRE_THROWS_AS(sstruct::Residue(""), sstruct::SecondaryStructureException);
        }
    }

}