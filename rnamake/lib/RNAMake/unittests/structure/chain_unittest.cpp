
#include "../common.hpp"

#include "base/file_io.h"
#include "base/settings.h"
#include "structure/chain.h"
#include "structure/is_equal.hpp"

TEST_CASE( "Test Chain for Structure", "[Chain]" ) {
    
    auto rts = ResidueTypeSet();
    auto path = unittest_resource_dir() + "structure/test_str_to_chain.dat";
    auto lines = get_lines_from_file(path);
    auto c = std::make_shared<Chain>(lines[0], rts);

    SECTION("can stringify chain") {
        auto s = c->to_str();
        auto c2 = std::make_shared<Chain>(s, rts);
        
        REQUIRE(are_chains_equal(c, c2, 0));
    }

    SECTION("can copy chain") {
        auto c2 = std::make_shared<Chain>(*c);
        REQUIRE(are_chains_equal(c, c2));
    }

    SECTION("testing sectioning chain into subchains using indices") {
        
        SECTION("cannot supply negative numbers") {
            REQUIRE_THROWS_AS(c->subchain(-1, 5), ChainException);
        }
        
        SECTION("end number is larger then the chain size") {
            REQUIRE_THROWS_AS(c->subchain(1, 10000), ChainException);
        }
        
        SECTION("subchain would be of size 0") {
            REQUIRE_THROWS_AS(c->subchain(1, 1), ChainException);
        }
        
        REQUIRE(c->subchain(0, 5)->length() == 5);
    }

    SECTION("testing sectioning chains using residues instead of indices") {
        
        auto r1 = c->get_residue(0);
        auto r2 = c->get_residue(5);
        
        REQUIRE(c->subchain(r1, r2)->length() == 5);
        
        auto r_copy = std::make_shared<Residue>(*r1);
        
        SECTION("using a residue not in the chain") {
            REQUIRE_THROWS_AS(c->subchain(r_copy, r2), ChainException);
        }
        
    }
    
    SECTION("stop empty chain from behaving weirdly") {
        auto c2 = std::make_shared<Chain>(ResidueOPs());
        
        REQUIRE_THROWS_AS(c2->first(), ChainException);
        REQUIRE_THROWS_AS(c2->last(), ChainException);

    }

    SECTION("get state") {
        auto cs = c->get_state();
        REQUIRE(cs->length() == c->length());

        auto cs_copy = std::make_shared<state::Chain>(*cs);
        REQUIRE(cs->length() == cs_copy->length());

        auto s = cs->to_str();
        cs_copy = std::make_shared<state::Chain>(s);
        REQUIRE(cs->length() == cs_copy->length());

    }
}









