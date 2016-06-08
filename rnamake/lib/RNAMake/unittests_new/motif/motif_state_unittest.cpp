
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "math/numerical.h"
#include "motif/motif_state.h"
#include "motif/motif.h"

TEST_CASE( "Test Motif states, motifs that dont have coordinates", "[MotifState]" ) {
    auto path = motif_dirs() + "base.motif";
    auto m = file_to_motif(path);
    path = motif_dirs() + "ref.motif";
    auto ref_m = file_to_motif(path);;

    auto ms = m->get_state();
    
    REQUIRE(ms->name() == m->name());
    REQUIRE(ms->score() == m->score());
    REQUIRE(ms->size() == m->residues().size());
    REQUIRE(ms->end_ids() == m->end_ids());
    REQUIRE(ms->block_end_add() == m->block_end_add());
    
    int i = -1 ;
    for(auto const & end : m->ends()) {
        i++;
        REQUIRE(end->name() == ms->end_names()[i]);
        REQUIRE(end->d() == ms->end_states()[i]->d());
        REQUIRE(are_xyzMatrix_equal(end->r(), ms->end_states()[i]->r()));
    }
    
    SECTION("test copy constructor") {
        auto ms_copy = std::make_shared<MotifState>(*ms);
        
        int i = -1 ;
        for(auto const & end : ms->end_states()) {
            i++;
            REQUIRE(ms_copy->end_names()[i] == ms->end_names()[i]);
            REQUIRE(ms_copy->end_states()[i]->d() == ms->end_states()[i]->d());
            REQUIRE(are_xyzMatrix_equal(ms_copy->end_states()[i]->r(),
                                        ms->end_states()[i]->r()));
        }
        
    }
    
    SECTION("test stringifying motif state") {
        auto s = ms->to_str();
        auto ms_copy = std::make_shared<MotifState>(s);
        
        int i = -1 ;
        for(auto const & end : ms->end_states()) {
            i++;
            REQUIRE(ms_copy->end_names()[i] == ms->end_names()[i]);
            REQUIRE(are_xyzVector_equal(ms_copy->end_states()[i]->d(),
                                        ms->end_states()[i]->d()));
            
            REQUIRE(are_xyzMatrix_equal(ms_copy->end_states()[i]->r(),
                                        ms->end_states()[i]->r()));
        }

        
        
    }
    
    SECTION("compare motif state align to motif align") {
        auto m2 = std::make_shared<Motif>(*m);
        auto m_aligned = get_aligned_motif(m->ends()[1], m2->ends()[0], m2);
        
        auto ms2 = std::make_shared<MotifState>(*ms);
        auto ms3 = std::make_shared<MotifState>(*ms);
        get_aligned_motif_state(ms->end_states()[1], ms2, ms3);
        

        REQUIRE(are_xyzVector_equal(m_aligned->ends()[1]->d(),
                                    ms2->end_states()[1]->d()));
        
        REQUIRE(are_xyzMatrix_equal(m_aligned->ends()[1]->r(),
                                    ms2->end_states()[1]->r()));
        
    }
    

}