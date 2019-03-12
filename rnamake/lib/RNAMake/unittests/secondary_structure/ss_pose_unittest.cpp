
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "secondary_structure/pose.h"
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Poses for secondary structure", "[SSPose]" ) {

    auto parser = secondary_structure::Parser();
    auto p = parser.parse_to_pose("GGGAGAAAACCCC",
                                  "(((.(....))))");
    
    REQUIRE(p->motifs().size() == 4);
    REQUIRE(p->sequence() == "GGGAGAAAACCCC");
    REQUIRE(p->dot_bracket() == "(((.(....))))");
    REQUIRE(p->ends().size() == 1);
    
    SECTION("test copy constructor for pose") {
        auto p_copy = std::make_shared<secondary_structure::Pose>(*p);
        
        REQUIRE(p_copy->ends().size() == 1);
        REQUIRE(p_copy->motifs().size() == 4);
        REQUIRE(p_copy->sequence() == "GGGAGAAAACCCC");
        REQUIRE(p_copy->dot_bracket() == "(((.(....))))");
        REQUIRE(p_copy->ends().size() == 1);
        
    }
    
    SECTION("test build helices") {
        auto helices = p->helices();
        REQUIRE(helices.size() == 1);
        REQUIRE(helices[0]->residues().size() == 6);
        
        parser.reset();
        auto p2 = parser.parse_to_pose("GGGAGGG+CCCCCC",
                                       "(((.(((+))))))");
        helices = p2->helices();
        REQUIRE(helices.size() == 2);
    }
    
    SECTION("test replacing the sequence of pose") {
        
        REQUIRE_THROWS_AS(p->replace_sequence("A"), secondary_structure::Exception);
        
        REQUIRE_THROWS_AS(p->replace_sequence("GGGAGAA&AACCCC"),
                          secondary_structure::Exception);
        
        
        p->replace_sequence("CCCAGAAAACGGG");
        REQUIRE(p->sequence() == "CCCAGAAAACGGG");
        
    }


}