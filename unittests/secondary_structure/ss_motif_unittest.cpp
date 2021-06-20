

//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "secondary_structure/motif.h"
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Motifs for secondary structure" ) {
    auto p = secondary_structure::Parser();
    auto m = p.parse_to_motif("GG+CC", "((+))");
    
    SUBCASE("test copy constructor") {
        auto m_copy = std::make_shared<secondary_structure::Motif>(*m);
        
        SUBCASE("all the residues have been copied correctly") {
            for(auto const & r : m->residues()) {
                CHECK(m_copy->get_residue(r->uuid()) != nullptr);
            }
        }
        
        SUBCASE("all the basepairs have been copied correctly") {
            for(auto const & bp : m->basepairs()) {
                CHECK(m_copy->get_basepair(bp->uuid()).size() != 0);
            }
        }
    }
    
    SUBCASE("test stringifing motif") {
        auto s = m->to_str();
        auto m_copy = std::make_shared<secondary_structure::Motif>(s);
        
        CHECK(m_copy->residues().size() == 4);
        CHECK(m_copy->basepairs().size() == 2);
        CHECK(m_copy->ends().size() == 2);
        CHECK(m_copy->end_ids().size() == 2);
    }

}