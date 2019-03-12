
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "motif/motif.h"
#include "resources/resource_manager.h"

TEST_CASE( "Test Resource Manager ", "[ResourceManager]" ) {
    
    SECTION("test construction of singleton") {
        REQUIRE_NOTHROW(RM::instance());
    }

    SECTION("test storing instance") {
        auto & rm = RM::instance();
        auto m = rm.motif("HELIX.IDEAL");
    }
    
    SECTION("test individual queries") {
        REQUIRE_NOTHROW(RM::instance().motif("HELIX.IDEAL"));
        REQUIRE_NOTHROW(RM::instance().motif("", "CC_LL_GG_RR"));
        REQUIRE_NOTHROW(RM::instance().motif("", "", "A5-B7"));
        REQUIRE_NOTHROW(RM::instance().motif("HELIX.IDEAL", "CC_LL_GG_RR"));
        REQUIRE_NOTHROW(RM::instance().motif("HELIX.IDEAL", "", "A5-B7"));
        REQUIRE_NOTHROW(RM::instance().motif("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7"));
        
        REQUIRE_NOTHROW(RM::instance().motif_state("HELIX.IDEAL"));
        REQUIRE_NOTHROW(RM::instance().motif_state("", "CC_LL_GG_RR"));
        REQUIRE_NOTHROW(RM::instance().motif_state("", "", "A5-B7"));
        REQUIRE_NOTHROW(RM::instance().motif_state("HELIX.IDEAL", "CC_LL_GG_RR"));
        REQUIRE_NOTHROW(RM::instance().motif_state("HELIX.IDEAL", "", "A5-B7"));
        REQUIRE_NOTHROW(RM::instance().motif_state("HELIX.IDEAL", "CC_LL_GG_RR", "A5-B7"));
        
        REQUIRE_NOTHROW(RM::instance().motif_state_ensemble("CC_LL_GG_RR"));

        REQUIRE_THROWS_AS(RM::instance().motif("TEST"), ResourceManagerException);
        REQUIRE_THROWS_AS(RM::instance().motif("", "TEST"), ResourceManagerException);
        REQUIRE_THROWS_AS(RM::instance().motif("", "", "TEST"), ResourceManagerException);

    }
    
    SECTION("test adding motif_state_ensembles") {
        
        auto path = base::unittest_resource_dir() + "resources/test.dat";
        REQUIRE_NOTHROW(RM::instance().register_extra_motif_ensembles(path));
        
        // should exist
        REQUIRE(RM::instance().has_supplied_motif_ensemble("TWOWAY.3TD0.1", "A7-B42") == 1);
        
        // does not exist wrong end_name
        REQUIRE(RM::instance().has_supplied_motif_ensemble("TWOWAY.3TD0.1", "A8-B42") == 0);

        auto me = RM::instance().get_supplied_motif_ensemble("TWOWAY.3TD0.1", "A7-B42");
        REQUIRE(me->size() == 1);
        
    }

    SECTION("test adding new motifs to resource manager") {
        
        REQUIRE_NOTHROW(RM::instance().add_motif(base::unittest_resource_dir() + "motif/GAAA_tetraloop"));

        auto m1 = RM::instance().motif("GAAA_tetraloop");
        REQUIRE(m1->name() == "GAAA_tetraloop");
        
    }
    
    
    
    
}