
//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

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
    
    SECTION("test pretty printing tree") {
        auto mt2 = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");

        mt2.add_motif(m1);
        mt2.add_motif(m2);
        //mt2.add_motif(m3, 1);
        //mt2.add_motif(nway);

        auto s = mt2.to_pretty_str();
        
        auto path = base::unittest_resource_dir() + "motif_tree/pretty_str_1.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto spl = base::split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            REQUIRE(spl[i] == lines[i-1]);
        }
        
    }
    
    SECTION("test pretty printing tree with branching") {
        auto mt2 = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        
        mt2.add_motif(m1);
        mt2.add_motif(nway);
        mt2.add_motif(m2);
        mt2.add_motif(m3, 1);
        
        auto s = mt2.to_pretty_str();
        
        auto path = base::unittest_resource_dir() + "motif_tree/pretty_str_2.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto spl = base::split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            REQUIRE(spl[i] == lines[i-1]);
        }
        
    }
    
    SECTION("test adding motifs to tree") {
        auto mt = MotifTree();
        
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        REQUIRE_NOTHROW(mt.add_motif(m1));
        REQUIRE(mt.size() == 1);
        
        REQUIRE_NOTHROW(mt.add_motif(m2));
        REQUIRE(mt.size() == 2);
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, 10), MotifTreeException);
        //REQUIRE_THROWS_AS(mt.add_motif(m3, 0), MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_motif(m3, 0, 10), MotifTreeException);
        
        SECTION("make sure motifs with the same uuid are rejected") {
            REQUIRE_THROWS_AS(mt.add_motif(m1), MotifTreeException);
        }
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, "A4-A5"), MotifTreeException);
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, "FAKE_END"), MotifTreeException);
        
    }
    
    SECTION("make sure weird chain topologies are caught") {
        auto hairpin1 = RM::instance().motif("HAIRPIN.4P95.2");
        auto hairpin2 = RM::instance().motif("HAIRPIN.4P95.2");
        auto mt = MotifTree();
        mt.set_option_value("sterics", false);
        mt.add_motif(hairpin1);
        mt.add_motif(hairpin2);

        REQUIRE_THROWS_AS(mt.get_structure(), MotifTreeException);
        REQUIRE_THROWS_AS(mt.secondary_structure(), MotifTreeException);
        REQUIRE_THROWS_AS(mt.to_pdb(), MotifTreeException);

    }
    
    SECTION("test getting nodes") {
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.3");
        auto m4 = RM::instance().motif("HELIX.IDEAL.3");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);

        REQUIRE_NOTHROW(mt.get_node(0));
        REQUIRE_NOTHROW(mt.get_node(m1->id()));
        REQUIRE_NOTHROW(mt.get_node("HELIX.IDEAL.3"));

        REQUIRE_THROWS_AS(mt.get_node(10), MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node("HELIX.IDEAL.2"), MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node("FAKE"), MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node(m4->id()), MotifTreeException);
    }
    
    SECTION("test remove motifs from tree") {
        
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        mt.add_motif(m1);
        mt.add_motif(m2);
        REQUIRE_NOTHROW(mt.remove_node(1));
        REQUIRE_THROWS_AS(mt.remove_node(1), MotifTreeException);
    }
    
    SECTION("test remove node levels from tree") {
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");

        mt.add_motif(m1);
        mt.increase_level();
        mt.add_motif(m2);
        mt.add_motif(m3);
        mt.remove_node_level();
        
        REQUIRE(mt.size() == 1);
        
        mt.add_motif(m2);
        mt.add_motif(m3);
        
        mt.remove_node_level(0);
        
        REQUIRE(mt.size() == 0);

        
    }
    
    SECTION("test loading motif tree from topology str") {
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);
        
        auto str = mt.topology_to_str();
        auto mt2 = MotifTree(str);
        
        REQUIRE(mt2.size() == 3);
        auto s = mt2.get_structure();
        
        REQUIRE(s->chains().size() == 2);
        REQUIRE(s->residues().size() > m1->residues().size());
        
    }
    
    SECTION("test connecting nodes") {
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        mt.add_motif(m1);
        mt.add_motif(nway);
        mt.add_motif(m2);
        
        mt.add_connection(1, 2, "", "");
        auto s = mt.get_structure();
        REQUIRE(s->chains().size() == 1);
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, 1), MotifTreeException);
        REQUIRE(mt.add_motif(m3) == -1);
        
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", ""), MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", m1->ends()[0]->name()), MotifTreeException);

    }
    
    SECTION("test building random motif trees") {
        auto g = helix_and_two_way();
        REQUIRE(g->size() == 2);
        
        auto builder = MotifTreeBuilder();
        auto mt = builder.build();
        REQUIRE(mt->size() == 2);
        REQUIRE(mt->get_node(0)->data()->mtype() == util::MotifType::HELIX);
        REQUIRE(mt->get_node(1)->data()->mtype() == util::MotifType::TWOWAY);
        
        auto mt2 = builder.build(3);
        REQUIRE(mt2->size() == 6);
        REQUIRE(mt2->get_node(2)->data()->mtype() == util::MotifType::HELIX);

    }
    
    SECTION("test building random motif tree capped with hairpin") {
        
        auto builder = MotifTreeBuilder(helix_and_two_way_and_hairpin());
        auto mt = builder.build();
        
        REQUIRE(mt->size() == 4);
        
    }
    
    SECTION("test more complex builds of random motif trees") {
        auto g = std::make_shared<BuilderGraph>();
        g->add_node(util::MotifType::HELIX, 2);
        g->add_node(util::MotifType::NWAY, 3);
        g->add_node(util::MotifType::HELIX, 2, 1, 1);
        g->add_node(util::MotifType::HELIX, 2, 1, 2);
        
        auto builder = MotifTreeBuilder(g);
        auto mt = builder.build(2);
        REQUIRE(mt->size() == 8);
    
        
    }
    
    SECTION("test copying") {
        auto mt2 = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        
        mt2.add_motif(m1);
        mt2.add_motif(nway);
        mt2.add_motif(m2);
        mt2.add_motif(m3, 1);

        mt2.add_connection(2, 3, "", "");
        
        auto mt_copy = MotifTree(mt2);
        
        REQUIRE(mt_copy.size() == mt2.size());
        
        auto m4 = RM::instance().motif("HELIX.IDEAL.2");
        REQUIRE(mt_copy.add_motif(m4) == -1);
        
        auto rna_struct = mt_copy.get_structure();
        REQUIRE(rna_struct->chains().size() == 1);
    }
    
}






















