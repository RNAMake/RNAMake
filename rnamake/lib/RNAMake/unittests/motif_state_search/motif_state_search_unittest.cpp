
//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_state_search/motif_state_search.h"



TEST_CASE( "Test Searching Motif States", "[MotifStateSearch]" ) {
    
    
    SECTION("test simple search") {
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("HELIX.IDEAL.3");
        auto m2 = RM::instance().motif("TWOWAY.2PN4.4");
        auto m3 = RM::instance().motif("HELIX.IDEAL.3");
        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3);
        
        auto start = mt.get_node(0)->data()->ends()[0]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto search = MotifStateSearch();
        search.set_option_value("accept_score", 0.5f);
        search.set_option_value("max_node_level", 4);
        search.set_option_value("verbose", false);
        
        search.setup(start, end);
        auto sol = search.next();

        REQUIRE(sol != nullptr);
        REQUIRE(sol->score() < 0.1);

    }
    
    SECTION("test miniTTR") {
        RM::instance().add_motif(unittest_resource_dir() + "motif/GAAA_tetraloop");
        
        auto mt = MotifTree();
        auto m1 = RM::instance().motif("GAAA_tetraloop", "", "A229-A245");
        auto m2 = RM::instance().motif("HELIX.IDEAL.3");
        auto m3 = RM::instance().motif("HELIX.IDEAL.3");

        mt.add_motif(m1);
        mt.add_motif(m2);
        mt.add_motif(m3, 0);
        
        auto start = mt.get_node(1)->data()->ends()[1]->state();
        auto end = mt.get_node(2)->data()->ends()[1]->state();
        auto beads = mt.beads();
        auto centers = Points();
        for(auto const & b : beads) {
            if(b.btype() != BeadType::PHOS) {
                centers.push_back(b.center());
            }
        }
        
        auto search = MotifStateSearch();
        search.set_option_value("verbose", false);
        search.beads(centers);
        //search.set_option_value("max_node_level", 6);
        search.setup(start, end);
        auto sol = search.next();
        REQUIRE(sol->score() < 10);
        
        
    }
    
}
