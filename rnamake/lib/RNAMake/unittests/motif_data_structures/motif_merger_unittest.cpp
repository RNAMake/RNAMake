
//headers for testing
#include "../common.hpp"

#include <sstream>

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "resources/motif_sqlite_library.h"
#include "motif_data_structures/motif_merger.h"

TEST_CASE( "Test Mergering Motifs into single structure ", "[MotifMerger]" ) {
    auto m1 = RM::instance().motif("HELIX.IDEAL.3");
    auto m2 = RM::instance().motif("HELIX.IDEAL.3");
    auto m3 = RM::instance().motif("HELIX.IDEAL.3");
    auto hairpin = RM::instance().motif("HAIRPIN.4P95.2");
    auto hairpin_2 = RM::instance().motif("HAIRPIN.4P95.2");
    
    SECTION("test adding two motifs in isolation") {
    
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2);
        auto s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 4);
        
    }
    
    SECTION("test adding two motifs in connected") {
        
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        auto s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 2);
        
    }
    
    SECTION("test adding helix and hairpin together for one chain") {
        
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(hairpin, hairpin->ends()[0], m1, m1->ends()[1]);
        auto s = mm.get_structure();

        REQUIRE(s->chains().size() == 1);
        
        mm = MotifMerger();
        mm.add_motif(hairpin);
        mm.add_motif(m1, m1->ends()[1], hairpin, hairpin->ends()[0]);
        s = mm.get_structure();
        
        REQUIRE(s->chains().size() == 1);
        
    }
    
    SECTION("test adding hairpins together for one chain") {
        // catch std::cout warnings
        auto ss = std::stringstream();
        auto old_buf = std::cout.rdbuf(ss.rdbuf());

        auto mm = MotifMerger();
        mm.add_motif(hairpin);
        mm.add_motif(hairpin_2, hairpin_2->ends()[0], hairpin, hairpin->ends()[0]);
        REQUIRE_THROWS_AS(mm.get_structure(), MotifMergerException);
        // warrnings should of been produced
        REQUIRE(ss.str().size() > 0);
        
        std::cout.rdbuf(old_buf); //reset
    }

    SECTION("catch build errors") {
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        REQUIRE_THROWS_AS(mm.add_motif(m3, m3->ends()[0], m1, m1->ends()[1]), MotifMergerException);
    }

    SECTION("catch updating motif that does not exist") {
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.update_motif(m1);
        REQUIRE_THROWS_AS(mm.update_motif(m2), MotifMergerException);
        
        
    }
    
    SECTION("catching motif with same id twice") {
        auto mm = MotifMerger();
        mm.add_motif(m1);
        REQUIRE_THROWS_AS(mm.add_motif(m1), MotifMergerException);
    }
 
    SECTION("test simple") {
        auto m1 = RM::instance().motif("HELIX.IDEAL");
        auto m2 = RM::instance().motif("HELIX.IDEAL");
        
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = MotifMerger();
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
        auto m1 = RM::instance().motif("HELIX.IDEAL.1");
        auto m2 = RM::instance().motif("TWOWAY.1GID.12");
        
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = MotifMerger();
        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        
        auto rna_struc = mm.get_structure();
        REQUIRE(rna_struc->get_residue(m1->ends()[0]->res1()->uuid()) != nullptr);
        REQUIRE(rna_struc->get_residue(m1->ends()[1]->res1()->uuid()) == nullptr);
     
        m1 = RM::instance().motif("HELIX.IDEAL.1");
        m2 = RM::instance().motif("TWOWAY.1GID.12");
        mm = MotifMerger();
        
        align_motif(m2->ends()[1]->state(), m1->ends()[0], m1);
        mm.add_motif(m2);
        mm.add_motif(m1, m1->ends()[0], m2, m2->ends()[1]);
        rna_struc = mm.get_structure();
        REQUIRE(rna_struc->get_residue(m1->ends()[0]->res1()->uuid()) == nullptr);
        REQUIRE(rna_struc->get_residue(m2->ends()[1]->res1()->uuid()) != nullptr);
    }
    
    SECTION("test sequence indenity conflict") {
        auto m1 = RM::instance().bp_step("GG_LL_CC_RR");
        auto m2 = RM::instance().bp_step("AA_LL_UU_RR");
        align_motif(m1->ends()[1]->state(), m2->ends()[0], m2);
        auto mm = MotifMerger();
        
        // catch std::cout warnings
        auto ss = std::stringstream();
        auto old_buf = std::cout.rdbuf(ss.rdbuf());

        mm.add_motif(m1);
        mm.add_motif(m2, m2->ends()[0], m1, m1->ends()[1]);
        
        REQUIRE(ss.str().size() > 0);
        std::cout.rdbuf(old_buf); //reset
    }
    
}









