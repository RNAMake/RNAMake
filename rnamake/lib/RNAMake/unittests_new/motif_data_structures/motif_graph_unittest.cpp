
//headers for testing
#include "../common.hpp"
#include "../tools/motif_graph_builder.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "base/settings.h"
#include "base/backtrace.hpp"
#include "resources/resource_manager.h"
#include "secondary_structure/util.h"
#include "structure/is_equal.hpp"
#include "motif_data_structures/motif_graph.h"



TEST_CASE( "Test Assembling Motifs together in Graph ", "[MotifGraph]" ) {

    
    SECTION("test adding motifs") {
        auto mg = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");

        mg.add_motif(m1);
        mg.add_motif(m2);
        mg.add_motif(m3);

        REQUIRE(mg.size() == 3);
    }
    
    SECTION("test removing motifs") {
        auto mg = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        
        mg.add_motif(m1);
        mg.add_motif(m2);
        
        mg.remove_motif(1);
        REQUIRE(mg.size() == 1);
        
        mg.add_motif(m2);
        mg.add_motif(m3);
        mg.remove_level(0);
        REQUIRE(mg.size() == 0);
        
    }
    
    SECTION("test getting nodes") {
        auto mg = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        mg.add_motif(m1);
        mg.add_motif(m2);
        
        REQUIRE_NOTHROW(mg.get_node(0));
        REQUIRE_THROWS(mg.get_node(10));
        
    }
    
    SECTION("test stringifying the topology of a motif graph") {
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(5);
        
        REQUIRE(mg->size() == 10);
        auto s = mg->topology_to_str();
        
        auto mg2 = std::make_shared<MotifGraph>(s, MotifGraphStringType::TOP);
        
        REQUIRE(mg->size() == mg2->size());
        
        for(int i = 0; i < mg->size(); i++) {
            auto atoms1 = mg->get_node(i)->data()->atoms();
            auto atoms2 = mg2->get_node(i)->data()->atoms();
            REQUIRE(are_atom_vectors_equal(atoms1, atoms2));
            
        }
        
        auto struc = mg2->get_structure();
        REQUIRE(struc->chains().size() == 2);
        
    }
    
    SECTION("test stringifying the motif graph yeilds indentical motifs when reloaded") {
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(5);
        auto s = mg->to_str();
        
        auto mg2 = std::make_shared<MotifGraph>(s, MotifGraphStringType::MG);

        REQUIRE(mg->size() == mg2->size());
        
        for(int i = 0; i < mg->size(); i++) {
            auto atoms1 = mg->get_node(i)->data()->atoms();
            auto atoms2 = mg2->get_node(i)->data()->atoms();
            REQUIRE(are_atom_vectors_equal(atoms1, atoms2));
            
        }
        
        auto struc = mg2->get_structure();
        REQUIRE(struc->chains().size() == 2);
    }
    
    SECTION("test compatibility with python stringification and topology serialization") {
        auto path = base_dir() + "/rnamake/unittests/resources/motif_graph/";
        auto lines = get_lines_from_file(path + "test.mg");
        auto mg = MotifGraph(lines[0], MotifGraphStringType::MG);
        auto s = mg.get_structure();
        
        lines = get_lines_from_file(path + "base_mg.mg");
        mg = MotifGraph(lines[0], MotifGraphStringType::MG);
        s = mg.get_structure();
        
    }
    
    SECTION("test copying motif graph") {
        auto builder = MotifGraphBuilder();
        auto mg = builder.build(5);
        
        auto mg2 = std::make_shared<MotifGraph>(*mg);
        
        REQUIRE(mg->size() == mg2->size());
        
        for(int i = 0; i < mg->size(); i++) {
            auto atoms1 = mg->get_node(i)->data()->atoms();
            auto atoms2 = mg2->get_node(i)->data()->atoms();
            REQUIRE(are_atom_vectors_equal(atoms1, atoms2));
            
        }

        
        
    }
 
    SECTION("test replacing idealized helices") {
        auto mg = MotifGraph();
        auto m = RM::instance().motif("HELIX.IDEAL.6");
        mg.add_motif(m);
        mg.replace_ideal_helices();
        
        REQUIRE(mg.size() == 7);
    
        auto builder = MotifGraphBuilder();
        auto mg2 = builder.build(2);
        
        SECTION("make sure complex build is producing the right number of motifs") {
        
            auto expected_total = 0;
            for(auto const & n : *mg2) {
                if(n->data()->mtype() != MotifType::HELIX) {
                    expected_total += 1;
                }
                else {
                    expected_total += n->data()->basepairs().size() -1;
                }
            }
        
            mg2->replace_ideal_helices();
            REQUIRE(mg2->size() == expected_total);
        }
        
        auto s = mg2->get_structure();
    }
    
    SECTION("test replacing helices with new sequence") {
        auto mg = MotifGraph();
        auto m = RM::instance().motif("HELIX.IDEAL.6");
        mg.add_motif(m);
        mg.replace_ideal_helices();

        auto dss = mg.designable_secondary_structure();
        sstruct::fill_basepairs_in_ss(dss);
        
        REQUIRE_NOTHROW(mg.replace_helical_sequence(dss));
        
        auto builder = MotifGraphBuilder();
        auto mg2 = builder.build(2);
        mg2->replace_ideal_helices();

        dss = mg2->designable_secondary_structure();
        sstruct::fill_basepairs_in_ss(dss);

        REQUIRE_NOTHROW(mg2->replace_helical_sequence(dss));

    }
    
    SECTION("test get end for easy building ") {
        auto base_path = base_dir() + "/rnamake/lib/RNAMake/apps/mini_ttr/resources/";
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");

        RM::instance().add_motif(base_path+"GAAA_tetraloop", "ttr");
        
        auto mg = MotifGraph();
        auto ttr_m = RM::instance().motif("ttr", "", "A229-A245");
        mg.add_motif(m1);
        mg.add_motif(ttr_m);
        mg.add_motif(m2);
        
        REQUIRE(mg.size() == 3);
        
        REQUIRE_NOTHROW(mg.get_available_end(1, "A222-A251"));

        REQUIRE_THROWS_AS(mg.get_available_end(0), MotifGraphException);
        REQUIRE_THROWS_AS(mg.get_available_end(2, "A4-A5"), MotifGraphException);

        SECTION("try getting end by the name of the motif and end name") {
        
            REQUIRE_NOTHROW(mg.get_available_end("ttr", "A222-A251"));
            REQUIRE_THROWS_AS(mg.get_available_end("ttr", "FAKE_END"), MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("FAKE_MOTIF", "A222-A251"), MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("HELIX.IDEAL.2", "A1-A8"), MotifGraphException);
            REQUIRE_THROWS_AS(mg.get_available_end("ttr", "A149-A154"), MotifGraphException);

        }
        
    }
    
    SECTION("test connecting motifs") {
        auto mg = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        mg.add_motif(m1);
        mg.add_motif(nway);
        mg.add_motif(m2);

        mg.add_connection(1, 2, "", "");
        auto s = mg.get_structure();
        REQUIRE(s->chains().size() == 1);

        REQUIRE_THROWS_AS(mg.add_motif(m3, -1, 1), MotifGraphException);
        REQUIRE(mg.add_motif(m3) == -1);
        
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "", ""), MotifGraphException);
        REQUIRE_THROWS_AS(mg.add_connection(1, 2, "", m1->ends()[0]->name()), MotifGraphException);
        
    }
    
    SECTION("test adding motif tree to motif graph") {
        
        auto g = std::make_shared<BuilderGraph>();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        g->add_node(HELIX, 2);
        g->add_node(NWAY, 3);
        g->add_node(HELIX, 2, 1, 1);
        g->add_node(HELIX, 2, 1, 2);
        
        auto builder = MotifTreeBuilder(g);
        auto mt = builder.build();
        
        auto mg = MotifGraph();
        mg.add_motif(m1);
        REQUIRE_NOTHROW(mg.add_motif_tree(mt));
        REQUIRE(mg.size() == 5);
        REQUIRE(mg.get_node(2)->data()->mtype() == MotifType::NWAY);
        
        auto mg2 = MotifGraph();
        mg2.add_motif(m1);
        REQUIRE_NOTHROW(mg2.add_motif_tree(mt, 0, "A1-A8"));
        REQUIRE(mg2.size() == 5);

        
    }
    
    SECTION("test pretty printing graph") {
        auto mg2 = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        
        mg2.add_motif(m1);
        mg2.add_motif(m2);
       
        
        auto s = mg2.to_pretty_str();
        
        auto path = unittest_resource_dir() + "motif_tree/pretty_str_1.dat";
        auto lines = get_lines_from_file(path);
        
        auto spl = split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            REQUIRE(spl[i] == lines[i-1]);
        }
        
    }
    
    SECTION("test pretty printing tree with branching") {
        auto mt2 = MotifGraph();
        auto m1 = RM::instance().motif("HELIX.IDEAL.2");
        auto m2 = RM::instance().motif("HELIX.IDEAL.2");
        auto m3 = RM::instance().motif("HELIX.IDEAL.2");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        
        mt2.add_motif(m1);
        mt2.add_motif(nway);
        mt2.add_motif(m2);
        mt2.add_motif(m3, 1);
        
        auto s = mt2.to_pretty_str();
        //std::cout << s << std::endl;
        
        auto path = unittest_resource_dir() + "motif_tree/pretty_str_2.dat";
        auto lines = get_lines_from_file(path);
        
        auto spl = split_str_by_delimiter(s, "\n");
        for(int i = 1; i < spl.size(); i++) {
            REQUIRE(spl[i] == lines[i-1]);
        }
        
    }
}



















