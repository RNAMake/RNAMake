

//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_tree.h"

#include <util/find_pair.h>

motif_data_structure::MotifStateTreeOP
_get_sub_tree() {
    auto mst = std::make_shared<motif_data_structure::MotifStateTree>();
    auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
    auto m2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
    auto m3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
    auto nway = resources::Manager::instance().motif_state("NWAY.1GID.0");
    mst->add_state(m1);
    mst->add_state(nway);
    mst->add_state(m2);
    mst->add_state(m3, 1);
    return mst;
}


TEST_CASE( "Test Assembling MotifStates together" ) {
    
    SUBCASE("test add states") {
        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto mst = motif_data_structure::MotifStateTree();
        REQUIRE_NOTHROW(mst.add_state(ms1));
        
        // cannot add the same state twice with same uuid
        REQUIRE_THROWS_AS(mst.add_state(ms1), motif_data_structure::MotifStateTreeException);
        
        // can never use parent_end_index=0 for a tree as that is where that node
        // is already connected to another node
        REQUIRE_THROWS_AS(mst.add_state(ms2, 0, 0), motif_data_structure::MotifStateTreeException);
        
        // catches invalid parent_index
        REQUIRE_THROWS_AS(mst.add_state(ms2, 0, 2), motif_data_structure::MotifStateTreeException);
        
        // invalid parent_end_name, is the name of end 0
        REQUIRE_THROWS_AS(mst.add_state(ms2, 0, "A4-A5"), motif_data_structure::MotifStateTreeException);
        
        // invalid parent_end_name, cannot be found as an end in state
        REQUIRE_THROWS_AS(mst.add_state(ms2, 0, "FAKE"), motif_data_structure::MotifStateTreeException);
        
        REQUIRE_NOTHROW(mst.add_state(ms2, 0, "A1-A8"));
        REQUIRE_NOTHROW(mst.add_state(ms3, -1, 1));
        CHECK(mst.size() == 3);
        
    }
    
    SUBCASE("test add another motif state tree") {
        auto mst = motif_data_structure::MotifStateTree();
        auto mst_add = _get_sub_tree();
        auto mst_add2 = _get_sub_tree();

        REQUIRE_NOTHROW(mst.add_mst(mst_add));
        CHECK(mst.size() == 4);
        
        REQUIRE_NOTHROW(mst.add_mst(mst_add2));
        CHECK(mst.size() == 8);

        auto mst2 = motif_data_structure::MotifStateTree();
        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        mst2.add_state(ms1);
        
        REQUIRE_NOTHROW(mst2.add_mst(mst_add, 0, "A1-A8"));
        CHECK(mst2.size() == 5);
        
    }
    
    SUBCASE("test adding connections") {
        auto mst = std::make_shared<motif_data_structure::MotifStateTree>();
        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif_state("NWAY.1GID.0", "", "A138-A180");
        mst->add_state(ms1);
        mst->add_state(nway);
        mst->add_state(ms2);
        
        // try connecting through 0th end position
        REQUIRE_THROWS_AS(mst->add_connection(1, 2, "A138-A180", ""), motif_data_structure::MotifStateTreeException);
        
        // try connecting thru an already used end position
        REQUIRE_THROWS_AS(mst->add_connection(1, 2, "A141-A162", ""), motif_data_structure::MotifStateTreeException);
        
        REQUIRE_NOTHROW(mst->add_connection(1, 2, "", ""));
        auto rna_struc = mst->get_structure();
        CHECK(rna_struc->chains().size() == 1);
        
        // cannot add new state where there is an connection present
        REQUIRE_THROWS_AS(mst->add_state(ms3, -1, 1), motif_data_structure::MotifStateTreeException);
        
        CHECK(mst->add_state(ms3) == -1);
        
    }
    
    SUBCASE("test copy") {
        auto mst2 = motif_data_structure::MotifStateTree();
        auto m1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto m3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto nway = resources::Manager::instance().motif_state("NWAY.1GID.0");
        
        mst2.add_state(m1);
        mst2.add_state(nway);
        mst2.add_state(m2);
        mst2.add_state(m3, 1);
        
        mst2.add_connection(2, 3, "", "");
        
        auto mst_copy = motif_data_structure::MotifStateTree(mst2);
        
        CHECK(mst_copy.size() == mst2.size());
        
        auto m4 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        CHECK(mst_copy.add_state(m4) == -1);
        
        auto rna_struct = mst_copy.get_structure();
        CHECK(rna_struct->chains().size() == 1);

    }
    
    SUBCASE("test generating new motif state tree from topology") {
        auto builder = MotifTreeBuilder();
        auto mt = builder.build(3);
        auto mst = std::make_shared<motif_data_structure::MotifStateTree>(mt);
        auto s = mst->topology_to_str();
        
        auto mst_copy = std::make_shared<motif_data_structure::MotifStateTree>(s);
        
        CHECK(mst->size() == mst_copy->size());
        
        int i = mst->size()-1;
        auto d1 = mst->get_node(i)->data()->get_end_state(1)->d();
        auto d2 = mst_copy->get_node(i)->data()->get_end_state(1)->d();
        
        CHECK(d1.distance(d2) < 0.1);

    }
    
    SUBCASE("test construction from motif_tree") {
        auto builder = MotifTreeBuilder();
        auto mt = builder.build(3);
        
        auto mst = std::make_shared<motif_data_structure::MotifStateTree>(mt);
        
        CHECK(mt->size() == mst->size());
        
        int i = mt->size()-1;
        auto d1 = mst->get_node(i)->data()->get_end_state(1)->d();
        auto d2 = mt->get_node(i)->data()->ends()[1]->d();
        
        CHECK(d1.distance(d2) < 0.1);
    }
    
    SUBCASE("test replacing a state in the tree") {
        //auto mt = builder.build(3);
        auto mst = std::make_shared<motif_data_structure::MotifStateTree>();
        mst->add_state(resources::Manager::instance().motif_state("HELIX.IDEAL.7"));
        mst->add_state(resources::Manager::instance().motif_state("HELIX.IDEAL.7"));
        mst->add_state(resources::Manager::instance().motif_state("HELIX.IDEAL.7"));

        auto new_state = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        mst->replace_state(2, new_state);

        auto mt2 = std::make_shared<motif_data_structure::MotifTree>();
        for(auto const & n : *mst) {            
            auto m = resources::Manager::instance().motif(n->data()->name(), "", n->data()->end_name(0));
            mt2->add_motif(m);
        }
        
        int i = mst->size()-1;
        auto d1 = mst->get_node(i)->data()->get_end_state(1)->d();
        auto d2 = mt2->get_node(i)->data()->ends()[1]->d();
        
        CHECK(d1.distance(d2) < 0.1);
    }

    SUBCASE("test removing nodes from tree") {
        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms2 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto ms3 = resources::Manager::instance().motif_state("HELIX.IDEAL.2");
        auto mst = motif_data_structure::MotifStateTree();
        mst.add_state(ms1);
        mst.add_state(ms2);
        mst.add_state(ms3);
        
        mst.remove_node();
        
        CHECK(mst.size() == 2);
        mst.add_state(ms3);

        mst.remove_node_level();
        CHECK(mst.size() == 0);
    
        mst.add_state(ms1);
        mst.increase_level();
        mst.add_state(ms2);
        mst.add_state(ms3);
        mst.remove_node_level();
        CHECK(mst.size() == 1);
        
    }
    
}



































