
//headers for testing
#include "../common.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_tree.h"

TEST_CASE( "Test Assembling Motifs together in Tree ", "[MotifTree]" ) {
    
    SECTION("test setting options") {
        auto mt = MotifTree();
        REQUIRE(mt.get_bool_option("sterics") == true);
        mt.set_option_value("sterics", false);
        REQUIRE(mt.get_bool_option("sterics") == false);
    }
    
    SECTION("test adding motifs to tree") {
        auto mt = MotifTree();
        
        auto m = RM::instance().motif("HELIX.IDEAL.2");
        REQUIRE_NOTHROW(mt.add_motif(m));
        REQUIRE(mt.size() == 1);
        
        REQUIRE_NOTHROW(mt.add_motif(m));
        REQUIRE(mt.size() == 2);
        
        REQUIRE_THROWS_AS(mt.add_motif(m, 10), MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_motif(m, 0), MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_motif(m, 0, 10), MotifTreeException);
        
        
        REQUIRE_THROWS_AS(mt.add_motif(m, -1, "FAKE_END"), MotifTreeException);
        
    }
    
    SECTION("make sure weird chain topologies are caught") {
        auto hairpin = RM::instance().motif("HAIRPIN.4P95.2");
        auto mt = MotifTree();
        mt.set_option_value("sterics", false);
        mt.add_motif(hairpin);
        mt.add_motif(hairpin);

        REQUIRE_THROWS_AS(mt.get_structure(), MotifTreeException);
        REQUIRE_THROWS_AS(mt.secondary_structure(), MotifTreeException);
        REQUIRE_THROWS_AS(mt.to_pdb(), MotifTreeException);

    }
    
    SECTION("test getting nodes") {
        auto mt = MotifTree();
        REQUIRE_THROWS_AS(mt.get_node(10), MotifTreeException);
    }
    
    SECTION("test remove motifs from tree") {
        
        auto mt = MotifTree();
        auto m = RM::instance().motif("HELIX.IDEAL.2");

        mt.add_motif(m);
        mt.add_motif(m);
        REQUIRE_NOTHROW(mt.remove_node(1));
        REQUIRE_THROWS_AS(mt.remove_node(1), MotifTreeException);
    }
    
    SECTION("test loading motif tree from topology str") {
        auto mt = MotifTree();
        auto m = RM::instance().motif("HELIX.IDEAL.2");
        mt.add_motif(m);
        mt.add_motif(m);
        mt.add_motif(m);
        
        auto str = mt.topology_to_str();
        auto mt2 = MotifTree(str);
        
        REQUIRE(mt2.size() == 3);
        auto s = mt2.get_structure();
        
        REQUIRE(s->chains().size() == 2);
        REQUIRE(s->residues().size() > m->residues().size());
        
    }
    
    SECTION("test connecting nodes") {
        auto mt = MotifTree();
        auto m = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        mt.add_motif(m);
        mt.add_motif(nway);
        mt.add_motif(m);
        
        mt.add_connection(1, 2, "", "");
        auto s = mt.get_structure();
        //REQUIRE(s->chains().size() == 1);
        
        REQUIRE_THROWS_AS(mt.add_motif(m, -1, 1), MotifTreeException);
        REQUIRE(mt.add_motif(m) == -1);
        
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", ""), MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", m->ends()[0]->name()), MotifTreeException);

    }
    

}