
//headers for testing
#include "../common.hpp"

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

    m2->new_res_uuids();
    m3->new_res_uuids();
    hairpin_2->new_res_uuids();
    
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
        
        auto mm = MotifMerger();
        mm.add_motif(hairpin);
        mm.add_motif(hairpin_2, hairpin_2->ends()[0], hairpin, hairpin->ends()[0]);
        REQUIRE_THROWS_AS(mm.get_structure(), MotifMergerException);
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
    
}