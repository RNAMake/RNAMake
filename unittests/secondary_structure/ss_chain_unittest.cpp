

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "util/uuid.h"
#include "secondary_structure/chain.h"

TEST_CASE( "Test Chain for Secondary Structure" ) {
    auto res =  secondary_structure::ResidueOPs();
    res.push_back(std::make_shared<secondary_structure::Residue>("G", "(", 1, "A", util::Uuid()));
    res.push_back(std::make_shared<secondary_structure::Residue>("A", ".", 2, "A", util::Uuid()));
    res.push_back(std::make_shared<secondary_structure::Residue>("C", ")", 3, "A", util::Uuid()));
    auto c = std::make_shared<secondary_structure::Chain>(res);
    
    CHECK(c->length() == 3);
    CHECK(c->sequence() == "GAC");
    CHECK(c->dot_bracket() == "(.)");
    
    SUBCASE("test copy constructor") {
        auto c2 = std::make_shared<secondary_structure::Chain>(*c);
        
        CHECK(c2->length() == 3);
        CHECK(c2->sequence() == "GAC");
        CHECK(c2->dot_bracket() == "(.)");
        
    }
    
    SUBCASE("test stringify of chain") {
        auto s = c->to_str();
        auto c2 = std::make_shared<secondary_structure::Chain>(s);

        CHECK(c2->length() == 3);
        CHECK(c2->sequence() == "GAC");
        CHECK(c2->dot_bracket() == "(.)");
    }
    
    SUBCASE("test errors from calling first and last on empty chain") {
        auto c = std::make_shared<secondary_structure::Chain>();
        
        REQUIRE_THROWS_AS(c->first(), secondary_structure::Exception);
        REQUIRE_THROWS_AS(c->last(), secondary_structure::Exception);
    }
    
}
