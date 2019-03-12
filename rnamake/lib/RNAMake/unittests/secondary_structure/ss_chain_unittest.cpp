


//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/chain.h"

TEST_CASE( "Test Chain for Secondary Structure", "[SSChain]" ) {
    auto res =  secondary_structure::ResidueOPs();
    res.push_back(std::make_shared<secondary_structure::Residue>("G", "(", 1, "A", util::Uuid()));
    res.push_back(std::make_shared<secondary_structure::Residue>("A", ".", 2, "A", util::Uuid()));
    res.push_back(std::make_shared<secondary_structure::Residue>("C", ")", 3, "A", util::Uuid()));
    auto c = std::make_shared<secondary_structure::Chain>(res);
    
    REQUIRE(c->length() == 3);
    REQUIRE(c->sequence() == "GAC");
    REQUIRE(c->dot_bracket() == "(.)");
    
    SECTION("test copy constructor") {
        auto c2 = std::make_shared<secondary_structure::Chain>(*c);
        
        REQUIRE(c2->length() == 3);
        REQUIRE(c2->sequence() == "GAC");
        REQUIRE(c2->dot_bracket() == "(.)");
        
    }
    
    SECTION("test stringify of chain") {
        auto s = c->to_str();
        auto c2 = std::make_shared<secondary_structure::Chain>(s);

        REQUIRE(c2->length() == 3);
        REQUIRE(c2->sequence() == "GAC");
        REQUIRE(c2->dot_bracket() == "(.)");
    }
    
    SECTION("test errors from calling first and last on empty chain") {
        auto c = std::make_shared<secondary_structure::Chain>();
        
        REQUIRE_THROWS_AS(c->first(), secondary_structure::Exception);
        REQUIRE_THROWS_AS(c->last(), secondary_structure::Exception);
    }
    
}
