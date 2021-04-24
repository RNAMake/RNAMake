
//headers for testing
#include "../common.hpp"

#include <sstream>

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "resources/motif_sqlite_library.h"
#include "motif_data_structure/motif_merger.h"

#include <util/find_pair.h>

TEST_CASE( "Test Mergering Motifs into single structure ", "[motif_data_structure::MotifMerger]" ) {
    auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.3");
    auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.3");
    auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.3");
    auto hairpin = resources::Manager::instance().motif("HAIRPIN.4P95.2");
    auto hairpin_2 = resources::Manager::instance().motif("HAIRPIN.4P95.2");
    
    SECTION("test adding two motifs in isolation") {
    
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2);
        auto s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 4);
        
    }
    
    SECTION("test adding two motifs in connected") {
        
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        auto s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 2);
        
    }
    
    SECTION("test adding helix and hairpin together for one chain") {
        
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(hairpin, hairpin->ends()[0], m1, m1->ends()[1]);
        auto s = mm.get_structure();

        REQUIRE(s->chains().size() == 1);
        
        mm = motif_data_structure::MotifMerger();
        mm.add_motif(hairpin);
        mm.add_motif(m1, m1->ends()[1], hairpin, hairpin->ends()[0]);
        s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 1);
        
    }
    
    SECTION("test adding hairpins together for one chain") {
        // catch std::cout warnings
        // We hit a most vexing parse if we use
        //auto ss(std::stringstream());
        // and use the deleted copy ctor if we use
        //auto ss = std::stringstream()
        //auto ss(std::stringstream(""));
        //std::stringstream ss;
        //auto old_buf = std::cout.rdbuf(ss.rdbuf());

        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(hairpin);
        mm.add_motif(hairpin_2, hairpin_2->ends()[0], hairpin, hairpin->ends()[0]);
        REQUIRE_THROWS_AS(mm.get_structure(), motif_data_structure::MotifMergerException);
        // warrnings should of been produced
        //REQUIRE(ss.str().size() > 0);
        
        //std::cout.rdbuf(old_buf); //reset
    }

    SECTION("catch build errors") {
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        REQUIRE_THROWS_AS(mm.add_motif(m3, m3->ends()[0], m1, m1->ends()[1]), motif_data_structure::MotifMergerException);
    }

    SECTION("catch updating motif that does not exist") {
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.update_motif(m1);
        REQUIRE_THROWS_AS(mm.update_motif(m2), motif_data_structure::MotifMergerException);
    }
    
    SECTION("catching motif with same id twice") {
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        REQUIRE_THROWS_AS(mm.add_motif(m1), motif_data_structure::MotifMergerException);
    }
 
    SECTION("test simple") {
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL");
        
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        
        auto rna_struc = mm.get_structure();
        REQUIRE(rna_struc->get_residue(m1->ends()[0]->res1()->uuid()) != nullptr);
        REQUIRE(rna_struc->get_residue(m1->ends()[1]->res1()->uuid()) == nullptr);

        auto ss = mm.secondary_structure();
        REQUIRE(ss->sequence() == "CCC&GGG");
        REQUIRE(ss->dot_bracket() == "(((&)))");

    }
    
    SECTION("test conserve sequence indentity with twoway") {
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.1");
        auto m2 = resources::Manager::instance().motif("TWOWAY.3P59.1");
        
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = motif_data_structure::MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        
        auto rna_struc = mm.get_structure();
        REQUIRE(rna_struc->get_residue(m1->ends()[0]->res1()->uuid()) != nullptr);
        REQUIRE(rna_struc->get_residue(m1->ends()[1]->res1()->uuid()) == nullptr);
     
        m1 = resources::Manager::instance().motif("HELIX.IDEAL.1");
        m2 = resources::Manager::instance().motif("TWOWAY.3P59.1");
        mm = motif_data_structure::MotifMerger();
        
        align_motif(m2->ends()[1]->state(), m1->ends()[0], m1);
        mm.add_motif(m2);
        mm.add_motif(m1, m1->ends()[0], m2, m2->ends()[1]);
        rna_struc = mm.get_structure();
        REQUIRE(rna_struc->get_residue(m1->ends()[0]->res1()->uuid()) == nullptr);
        REQUIRE(rna_struc->get_residue(m2->ends()[1]->res1()->uuid()) != nullptr);
    }
    
    SECTION("test sequence indenity conflict") {
        auto m1 = resources::Manager::instance().bp_step("GG_LL_CC_RR");
        auto m2 = resources::Manager::instance().bp_step("AA_LL_UU_RR");
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = motif_data_structure::MotifMerger();
        
        // catch std::cout warnings
        //auto ss(std::stringstream(""));
        //std::stringstream ss;
        //auto old_buf = std::cout.rdbuf(ss.rdbuf());

        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        
        //REQUIRE(ss.str().size() > 0);
        //std::cout.rdbuf(old_buf); //reset
    }
    
}









