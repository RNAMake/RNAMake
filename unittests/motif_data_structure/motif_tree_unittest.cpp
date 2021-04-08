

//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_tree.h"

TEST_CASE( "Test Assembling Motifs together in Tree " ) {
    
    SUBCASE("test setting options") {
        auto mt = motif_data_structure::MotifTree();
        CHECK(mt.get_bool_option("sterics") == true);
        mt.set_option_value("sterics", false);
        CHECK(mt.get_bool_option("sterics") == false);
    }
    
    SUBCASE("test pretty printing tree") {
        auto mt2 = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");

        mt2.add_motif(m1);
        mt2.add_motif(m2);
        //mt2.add_motif(m3, 1);
        //mt2.add_motif(nway);

        auto s = mt2.to_pretty_str();
        
        auto path = base::unittest_resource_dir() + "motif_tree/pretty_str_1.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto spl = base::split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            CHECK(spl[i] == lines[i-1]);
        }
        
    }
    
    SUBCASE("test pretty printing tree with branching") {
        auto mt2 = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");
        
        mt2.add_motif(m1);
        mt2.add_motif(nway);
        mt2.add_motif(m2);
        mt2.add_motif(m3, 1);
        
        auto s = mt2.to_pretty_str();
        
        auto path = base::unittest_resource_dir() + "motif_tree/pretty_str_2.dat";
        auto lines =base::get_lines_from_file(path);
        
        auto spl = base::split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            CHECK(spl[i] == lines[i-1]);
        }
        
    }
    
    SUBCASE("test adding motifs to tree") {
        auto mt = motif_data_structure::MotifTree();
        
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        REQUIRE_NOTHROW(mt.add_motif(m1));
        CHECK(mt.size() == 1);
        
        REQUIRE_NOTHROW(mt.add_motif(m2));
        CHECK(mt.size() == 2);
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, 10), motif_data_structure::MotifTreeException);
        //REQUIRE_THROWS_AS(mt.add_motif(m3, 0), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_motif(m3, 0, 10), motif_data_structure::MotifTreeException);
        
        SUBCASE("make sure motifs with the same uuid are rejected") {
            REQUIRE_THROWS_AS(mt.add_motif(m1), motif_data_structure::MotifTreeException);
        }
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, "A4-A5"), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, "FAKE_END"), motif_data_structure::MotifTreeException);
        
    }
    
    SUBCASE("make sure weird chain topologies are caught") {
        auto hairpin1 = resources::Manager::instance().motif("HAIRPIN.4P95.2");
        auto hairpin2 = resources::Manager::instance().motif("HAIRPIN.4P95.2");
        auto mt = motif_data_structure::MotifTree();
        mt.set_option_value("sterics", false);
        mt.add_motif(hairpin1);
        mt.add_motif(hairpin2);

        REQUIRE_THROWS_AS(mt.get_structure(), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.secondary_structure(), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.to_pdb(), motif_data_structure::MotifTreeException);

    }
    
    SUBCASE("test getting nodes") {
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.3");
        auto m4 = resources::Manager::instance().motif("HELIX.IDEAL.3");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);

        REQUIRE_NOTHROW(mt.get_node(0));
        REQUIRE_NOTHROW(mt.get_node(m1->id()));
        REQUIRE_NOTHROW(mt.get_node("HELIX.IDEAL.3"));

        REQUIRE_THROWS_AS(mt.get_node(10), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node("HELIX.IDEAL.2"), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node("FAKE"), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.get_node(m4->id()), motif_data_structure::MotifTreeException);
    }
    
    SUBCASE("test remove motifs from tree") {
        
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        mt.add_motif(m1);
        mt.add_motif(m2);
        REQUIRE_NOTHROW(mt.remove_node(1));
        REQUIRE_THROWS_AS(mt.remove_node(1), motif_data_structure::MotifTreeException);
    }
    
    SUBCASE("test remove node levels from tree") {
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");

        mt.add_motif(m1);
        mt.increase_level();
        mt.add_motif(m2);
        mt.add_motif(m3);
        mt.remove_node_level();
        
        CHECK(mt.size() == 1);
        
        mt.add_motif(m2);
        mt.add_motif(m3);
        
        mt.remove_node_level(0);
        
        CHECK(mt.size() == 0);

        
    }
    
    SUBCASE("test loading motif tree from topology str") {
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);
        
        auto str = mt.topology_to_str();
        auto mt2 = motif_data_structure::MotifTree(str);
        
        CHECK(mt2.size() == 3);
        auto s = mt2.get_structure();
        
        CHECK(s->chains().size() == 2);
        CHECK(s->residues().size() > m1->residues().size());
        
    }
    
    SUBCASE("test connecting nodes") {
        auto mt = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");
        mt.add_motif(m1);
        mt.add_motif(nway);
        mt.add_motif(m2);
        
        mt.add_connection(1, 2, "", "");
        auto s = mt.get_structure();
        CHECK(s->chains().size() == 1);
        
        REQUIRE_THROWS_AS(mt.add_motif(m3, -1, 1), motif_data_structure::MotifTreeException);
        CHECK(mt.add_motif(m3) == -1);
        
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", ""), motif_data_structure::MotifTreeException);
        REQUIRE_THROWS_AS(mt.add_connection(1, 2, "", m1->ends()[0]->name()), motif_data_structure::MotifTreeException);

    }
    
    SUBCASE("test building random motif trees") {
        auto g = helix_and_two_way();
        CHECK(g->size() == 2);
        
        auto builder = MotifTreeBuilder();
        auto mt = builder.build();
        CHECK(mt->size() == 2);
        CHECK(mt->get_node(0)->data()->mtype() == util::MotifType::HELIX);
        CHECK(mt->get_node(1)->data()->mtype() == util::MotifType::TWOWAY);
        
        auto mt2 = builder.build(3);
        CHECK(mt2->size() == 6);
        CHECK(mt2->get_node(2)->data()->mtype() == util::MotifType::HELIX);

    }
    
    SUBCASE("test building random motif tree capped with hairpin") {
        
        auto builder = MotifTreeBuilder(helix_and_two_way_and_hairpin());
        auto mt = builder.build();
        
        CHECK(mt->size() == 4);
        
    }
    
    SUBCASE("test more complex builds of random motif trees") {
        auto g = std::make_shared<BuilderGraph>();
        g->add_node(util::MotifType::HELIX, 2);
        g->add_node(util::MotifType::NWAY, 3);
        g->add_node(util::MotifType::HELIX, 2, 1, 1);
        g->add_node(util::MotifType::HELIX, 2, 1, 2);
        
        auto builder = MotifTreeBuilder(g);
        auto mt = builder.build(2);
        CHECK(mt->size() == 8);
    
        
    }
    
    SUBCASE("test copying") {
        auto mt2 = motif_data_structure::MotifTree();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");
        
        mt2.add_motif(m1);
        mt2.add_motif(nway);
        mt2.add_motif(m2);
        mt2.add_motif(m3, 1);

        mt2.add_connection(2, 3, "", "");
        
        auto mt_copy = motif_data_structure::MotifTree(mt2);
        
        CHECK(mt_copy.size() == mt2.size());
        
        auto m4 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        CHECK(mt_copy.add_motif(m4) == -1);
        
        auto rna_struct = mt_copy.get_structure();
        CHECK(rna_struct->chains().size() == 1);
    }
    
}






















