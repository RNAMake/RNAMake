
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "secondary_structure/motif.h"
#include "secondary_structure/secondary_structure_parser.h"

TEST_CASE( "Test Motifs for secondary structure", "[SSMotif]" ) {
    auto p = sstruct::SecondaryStructureParser();
    auto m = p.parse_to_motif("GG+CC", "((+))");
    
    SECTION("test copy constructor") {
        auto m_copy = std::make_shared<sstruct::Motif>(*m);
        
        SECTION("all the residues have been copied correctly") {
            for(auto const & r : m->residues()) {
                REQUIRE(m_copy->get_residue(r->uuid()) != nullptr);
            }
        }
        
        SECTION("all the basepairs have been copied correctly") {
            for(auto const & bp : m->basepairs()) {
                REQUIRE(m_copy->get_basepair(bp->uuid()).size() != 0);
            }
        }
    }
    
    SECTION("test stringifing motif") {
        auto s = m->to_str();
        auto m_copy = std::make_shared<sstruct::Motif>(s);
        
        REQUIRE(m_copy->residues().size() == 4);
        REQUIRE(m_copy->basepairs().size() == 2);
        REQUIRE(m_copy->ends().size() == 2);
        REQUIRE(m_copy->end_ids().size() == 2);
    }

}