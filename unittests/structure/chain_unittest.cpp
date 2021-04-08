

#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/chain.h"
#include "structure/is_equal.h"

TEST_CASE( "Test Chain for Structure" ) {
    
    auto rts = structure::ResidueTypeSet();
    auto path = base::unittest_resource_dir() + "structure/test_str_to_chain.dat";
    auto lines =base::get_lines_from_file(path);
    auto c = std::make_shared<structure::Chain>(lines[0], rts);
    
    SUBCASE("can stringify chain") {
        auto s = c->to_str();
        auto c2 = std::make_shared<structure::Chain>(s, rts);
        
        CHECK(are_chains_equal(c, c2, 0));
    }
    
    SUBCASE("can copy chain") {
        auto c2 = std::make_shared<structure::Chain>(*c);
        
        CHECK(are_chains_equal(c, c2));
    }
    
    SUBCASE("testing sectioning chain into subchains using indices") {
        
        SUBCASE("cannot supply negative numbers") {
            REQUIRE_THROWS_AS(c->subchain(-1, 5), structure::ChainException);
        }
        
        SUBCASE("end number is larger then the chain size") {
            REQUIRE_THROWS_AS(c->subchain(1, 10000), structure::ChainException);
        }
        
        SUBCASE("subchain would be of size 0") {
            REQUIRE_THROWS_AS(c->subchain(1, 1), structure::ChainException);
        }
        
        CHECK(c->subchain(0, 5)->length() == 5);
    }
    
    SUBCASE("testing sectioning chains using residues instead of indices") {
        
        auto r1 = c->residues()[0];
        auto r2 = c->residues()[5];
        
        CHECK(c->subchain(r1, r2)->length() == 5);
        
        auto r_copy = std::make_shared<structure::Residue>(*r1);
        
        SUBCASE("using a residue not in the chain") {
            REQUIRE_THROWS_AS(c->subchain(r_copy, r2), structure::ChainException);
        }
        
    }
    
    SUBCASE("stop empty chain from behaving weirdly") {
        auto c2 = std::make_shared<structure::Chain>();
        
        REQUIRE_THROWS_AS(c2->first(), structure::ChainException);
        REQUIRE_THROWS_AS(c2->last(), structure::ChainException);

    }
}









