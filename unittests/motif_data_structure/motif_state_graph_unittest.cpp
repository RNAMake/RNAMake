//
//  motif_state_graph_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//


#include <stdio.h>

//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_graph.hpp"

TEST_CASE( "Test Assembling MotifStates together in a graph " ) {
    
    SUBCASE("test adding states") {
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        
        SUBCASE("cannot find end if there is not parent") {
            REQUIRE_THROWS_AS(msg->add_state(ms1, -1, "A1-A8"), motif_data_structure::MotifStateGraphException);
        }
        
        msg->add_state(ms1);
        
        CHECK(msg->size() == 1);
        
        // can never use parent_end_index=0 for a graph as that is where that node
        // is already connected to another node
        REQUIRE_THROWS_AS(msg->add_state(ms2, -1, 0), motif_data_structure::MotifStateGraphException);
        
        // catches invalid parent_index
        REQUIRE_THROWS_AS(msg->add_state(ms2, 2), motif_data_structure::MotifStateGraphException);
        
        // invalid parent_end_index, has only 0 and 1
        REQUIRE_THROWS_AS(msg->add_state(ms2, -1, 3), motif_data_structure::MotifStateGraphException);
        
        // invalid parent_end_name, is the name of end 0
        REQUIRE_THROWS_AS(msg->add_state(ms2, -1, "A4-A5"), motif_data_structure::MotifStateGraphException);
        
        // invalid parent_end_name, cannot be found as an end in states
        REQUIRE_THROWS_AS(msg->add_state(ms2, -1, "FAKE"), motif_data_structure::MotifStateGraphException);

    }
 
    SUBCASE("test setup from mg") {
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto tc = resources::Manager::instance().motif("TC.1S72.0");
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        mg->add_motif(tc);
        mg->add_motif(m1);
        mg->add_motif(m2, 0);
        mg->add_connection(0, 2, "", "");
        
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        CHECK(mg->size() == msg->size());
        
    }
    
    SUBCASE("test to motif graph") {
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto tc = resources::Manager::instance().motif("TC.1S72.0");
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        mg->add_motif(tc);
        mg->add_motif(m1);
        mg->add_motif(m2, 0);
        
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        auto mg2 = msg->to_motif_graph();
        
        CHECK(mg->size() == mg2->size());
        
        auto atoms1 = mg->get_structure()->atoms();
        auto atoms2 = mg2->get_structure()->atoms();
        
        for(int i = 0; i < atoms1.size(); i++) {
            auto diff = atoms1[i]->get_coords().distance(atoms2[i]->get_coords());
            CHECK(diff < 0.1);
        }
        
        mg->add_connection(0, 2, "", "");
        msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        mg2 = msg->to_motif_graph();
        
        CHECK(mg->get_structure()->residues().size() == mg2->get_structure()->residues().size());
        
        
    }
    
    SUBCASE("test to motif graph 2") {
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL.2");
        m2->move(math::Point{40, 0, 0});
        
        auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        mg->add_motif(m1);
        mg->add_motif(m2, -1, -1, 1);
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"));
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"));
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"));
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"), 0);
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"));
        mg->add_motif(resources::Manager::instance().motif("HELIX.IDEAL.2"));
        
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        auto mg2 = msg->to_motif_graph();
        CHECK(mg->size() == mg2->size());
        
        mg->add_connection(4, 7, "", "");
        msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        mg2 = msg->to_motif_graph();
        
        CHECK(mg2->get_structure()->chains().size() == 2);
    }
    
    SUBCASE("test adding connections") {
        auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto tc = resources::Manager::instance().motif("TC.1S72.0")->get_state();
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->add_state(tc);
        msg->add_state(m1);
        msg->add_state(m2, 0);
        
        CHECK(msg->size() == 3);
        
        msg->add_connection(0, 2, "", "");
        auto mg = msg->to_motif_graph();
        
        CHECK(mg->get_structure()->chains().size() == 1);

    }
    
    SUBCASE("test remove state") {
        auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->add_state(m1);
        msg->add_state(m2);
        msg->add_state(m3);
        
        msg->remove_state(2);
        CHECK(msg->size() == 2);
        
        msg->add_state(m3);
        
        auto mg = msg->to_motif_graph();
        CHECK(mg->get_structure()->chains().size() == 2);
        

    }
    
    SUBCASE("test remove node level") {
        auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg->add_state(m1);
        msg->increase_level();
        msg->add_state(m2);
        msg->add_state(m3);
        msg->remove_level(1);
        CHECK(msg->size() == 1);
    }
    
    SUBCASE("test replace state") {
        auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif_state("TWOWAY.2PN4.4");
        auto m3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        
        msg->add_state(m1);
        msg->add_state(m2);
        msg->add_state(m3);
        
        msg->replace_state(1, resources::Manager::instance().motif_state("HELIX.IDEAL.2"));
        
        CHECK(msg->get_node(1)->data()->name() == "HELIX.IDEAL.2");
        
        auto msg2 = std::make_shared<motif_data_structure::MotifStateGraph>();
        msg2->add_state(m1);
        msg2->add_state(m3);
        msg2->add_state(resources::Manager::instance().motif_state("HELIX.IDEAL.2"));
        
        auto d1 = msg->get_node(2)->data()->get_end_state(1)->d();
        auto d2 = msg2->get_node(2)->data()->get_end_state(1)->d();
        CHECK(d1.distance(d2) < 0.001);
        
    }
    
    /*SUBCASE("test multiple alignments") {
        auto path = base::base_dir() + "/rnamake/unittests/resources/motif_graph/";
        auto lines =base::get_lines_from_file(path+"tecto_chip_only.mg");
        auto mg = std::make_shared<motif_data_structure::MotifGraph>(lines[0], motif_data_structure::MotifGraphStringType::MG);
        auto msg = std::make_shared<motif_data_structure::MotifStateGraph>(mg);
        auto mg2 = msg->to_motif_graph();
        
        //CHECK(mg->sequence() == mg2->sequence());
        //std::cout << mg2->sequence() << std::endl;

    }*/
    
}
