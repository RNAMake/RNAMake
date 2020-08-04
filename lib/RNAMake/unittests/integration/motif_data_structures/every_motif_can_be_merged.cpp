
//headers for testing
#include "../../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "resources/motif_sqlite_library.h"
#include "motif_data_structure/motif_merger.h"

TEST_CASE( "Test Mergering Motifs into single structure ", "[motif_data_structure::MotifMerger]" ) {
    auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.3");
    auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.3");
    m2->new_res_uuids();
    
    SECTION("test merging all two way junctions with flanking helices") {
        
        auto mlib = resources::resources::MotifSqliteLibrary("twoway");
        mlib.load_all();
        
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        int res_count = 0;
        
        for(auto const & m : mlib) {
            mm.add_motif(m, m->ends()[0], m1, m1->ends()[1]);
            mm.add_motif(m2, m2->ends()[0], m, m->ends()[1]);
            res_count = (int)m->residues().size() + (int)m2->residues().size() +
            (int)m1->residues().size() - 4;
            
            auto s = mm.get_structure();
            auto ss = mm.secondary_structure();
            
            /*auto m1_ss = m1->secondary_structure();
             auto m_ss = m->secondary_structure();
             
             auto seq1 = m1_ss->chains()[0]->sequence().substr(0, m1_ss->chains()[0]->length()-1) +
             m_ss->chains()[0]->sequence();
             auto seq2 = m_ss->chains()[1]->sequence() + m1_ss->chains()[1]->sequence().substr(1);
             
             std::cout << seq1 << "&" << seq2 << std::endl;
             */
            
            //std::cout << m2->sequence() << " " << m2->residues().size() << std::endl;
            
            REQUIRE(s->chains().size() == 2);
            REQUIRE(s->residues().size() == res_count);
            
            mm.remove_motif(m);
            mm.remove_motif(m2);
            
        }
    }
    
    SECTION("test merging all hairpins with a helix") {
        auto mlib = resources::resources::MotifSqliteLibrary("hairpin");
        mlib.load_all();
        
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        int res_count = 0;

        for(auto const & m : mlib) {
            
            mm.add_motif(m, m->ends()[0], m1, m1->ends()[1]);
            auto s = mm.get_structure();

            REQUIRE(s->chains().size() == 1);

            mm.remove_motif(m);
        }
        
    }
    
    SECTION("test merging all nway junctions") {
        auto mlib = resources::resources::MotifSqliteLibrary("nway");
        mlib.load_all();
        
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        int res_count = 0;
        
        for(auto const & m : mlib) {
            
            mm.add_motif(m, m->ends()[0], m1, m1->ends()[1]);
            mm.add_motif(m2, m2->ends()[0], m, m->ends()[1]);
            auto s = mm.get_structure();
            
            //std::cout << m->name() <<  " " << m->ends()[0]->name() <<  std::endl;
            REQUIRE(s->chains().size() == m->chains().size());
            
            mm.remove_motif(m);
            mm.remove_motif(m2);
        }

        
        
    }
    
}

















