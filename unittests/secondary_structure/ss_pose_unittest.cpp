

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include <secondary_structure/pose.h>
#include <secondary_structure/secondary_structure_parser.h>
#include <secondary_structure/sequence_tools.h>
#include <secondary_structure/sequence_constraint.h>

TEST_CASE( "Test Poses for secondary structure" ) {

    auto parser = secondary_structure::Parser();
    auto p = parser.parse_to_pose("GGGAGAAAACCCC",
                                  "(((.(....))))");
    
    CHECK(p->motifs().size() == 4);
    CHECK(p->sequence() == "GGGAGAAAACCCC");
    CHECK(p->dot_bracket() == "(((.(....))))");
    CHECK(p->ends().size() == 1);
    
    SUBCASE("test copy constructor for pose") {
        auto p_copy = std::make_shared<secondary_structure::Pose>(*p);
        
        CHECK(p_copy->ends().size() == 1);
        CHECK(p_copy->motifs().size() == 4);
        CHECK(p_copy->sequence() == "GGGAGAAAACCCC");
        CHECK(p_copy->dot_bracket() == "(((.(....))))");
        CHECK(p_copy->ends().size() == 1);
        
    }
    
    SUBCASE("test build helices") {
        auto helices = p->helices();
        CHECK(helices.size() == 1);
        CHECK(helices[0]->residues().size() == 6);
        
        parser.reset();
        auto p2 = parser.parse_to_pose("GGGAGGG+CCCCCC",
                                       "(((.(((+))))))");
        helices = p2->helices();
        CHECK(helices.size() == 2);
    }
    
    SUBCASE("test replacing the sequence of pose") {
        
        REQUIRE_THROWS_AS(p->replace_sequence("A"), secondary_structure::Exception);
        
        REQUIRE_THROWS_AS(p->replace_sequence("GGGAGAA&AACCCC"),
                          secondary_structure::Exception);
        
        
        p->replace_sequence("CCCAGAAAACGGG");
        CHECK(p->sequence() == "CCCAGAAAACGGG");
    }

    SUBCASE("test sequence tools") {
        SUBCASE("test convert_res_name_to_type()") {
            CHECK(secondary_structure::convert_res_name_to_type('A') == secondary_structure::ResType::ADE);
            REQUIRE_THROWS_AS(secondary_structure::convert_res_name_to_type('F'), secondary_structure::Exception);
        }

        SUBCASE("test get_res_types_from_sequence()") {
            auto str = String("AUCG");
            auto res_codes = secondary_structure::ResTypes();
            secondary_structure::get_res_types_from_sequence(str, res_codes);
            CHECK(res_codes.size() == 4);
            CHECK(res_codes[0] == secondary_structure::ResType::ADE);
            CHECK(res_codes[1] == secondary_structure::ResType::URA);
            CHECK(res_codes[2] == secondary_structure::ResType::CYT);
            CHECK(res_codes[3] == secondary_structure::ResType::GUA);
        }

        SUBCASE("test find_res_types_in_pose()") {
            auto str = String("CCCC");
            auto res_codes = secondary_structure::ResTypes();
            secondary_structure::get_res_types_from_sequence(str, res_codes);
            CHECK(secondary_structure::find_res_types_in_pose(p, res_codes) == 1);

            auto p_new = parser.parse_to_pose("AAAAGAAAAUUUU",
                                              "(((.(....))))");

            // multiple instances
            str = String("AAA");
            res_codes = secondary_structure::ResTypes();
            secondary_structure::get_res_types_from_sequence(str, res_codes);
            CHECK(secondary_structure::find_res_types_in_pose(p_new, res_codes) == 2);

            // target array too long
            str = String("AAAAAAAAAAAAAAAAAAAA");
            res_codes = secondary_structure::ResTypes();
            secondary_structure::get_res_types_from_sequence(str, res_codes);
            CHECK(secondary_structure::find_res_types_in_pose(p_new, res_codes) == 0);
        }

        SUBCASE("test find_gc_helix_stretches()") {
            auto stretches = secondary_structure::find_gc_helix_stretches(p, 3);
            CHECK(stretches == 1);

            auto p_new = parser.parse_to_pose("GCGAUAUGGGUCGAGCCCAAGUUAGGGAAACCUAGAGGGCAGUGAAAGACCCUAAGUCGC",
                                              "(((((.((((((..((((....((((....))))..)))).......)))))..))))))");

            CHECK(secondary_structure::find_gc_helix_stretches(p_new, 3) == 3);
        }

    }

    SUBCASE("test sequence constraints") {
        auto seq_constraint = secondary_structure::DisallowedSequence("CCCC");
        CHECK(seq_constraint.violations(p) == 1);

        auto seq_constraint_2 = secondary_structure::GCHelixStretchLimit(3);
        CHECK(seq_constraint_2.violations(p) == 1);

        auto seq_constraints = secondary_structure::SequenceConstraints();
        seq_constraints.add_disallowed_sequence("CCCC");
        seq_constraints.add_gc_helix_stretch_limit(3);
        auto violations = seq_constraints.violations(p);
        CHECK(violations.size() == 2);
        CHECK(violations[0] == 1);
        CHECK(violations[1] == 1);

    }

}






















