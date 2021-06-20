//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_data_structure/motif_state_ensemble_tree.h"
#include "motif_data_structure/motif_state_ensemble_graph.h"

#include <util/find_pair.h>


TEST_CASE( "Test Assembling MotifEnsembleStates together" ) {
    
    SUBCASE("test adding ensembles to tree") {
    
        auto mset = motif_data_structure::MotifStateEnsembleTree();
        auto mse = resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR");
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        
        CHECK(mset.size() == 2);
    }
    
    SUBCASE("test getting a motif state tree from the ensemble tree") {
        
        auto mset = motif_data_structure::MotifStateEnsembleTree();
        auto mse = resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR");
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);

        auto mst = mset.to_mst();
        CHECK(mst->size() == mset.size());
        
    }
    
    SUBCASE("setup from motif tree") {
        auto mt = std::make_shared<motif_data_structure::MotifTree>();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL");
        mt->add_motif(m1);
        mt->add_motif(m2);

        auto ms1 = resources::Manager::instance().motif_state("HELIX.IDEAL");
        auto mset = motif_data_structure::MotifStateEnsembleTree(mt);
        
        CHECK(mset.size() == 2);
        
        auto mst = mset.to_mst();
        //CHECK(mst->get_node(0)->data()->name() == "CG=CG.0");
        
        auto mt2 = std::make_shared<motif_data_structure::MotifTree>();
        auto m3 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto nway = resources::Manager::instance().motif("NWAY.1GID.0");
        
        mt2->add_motif(m1);
        mt2->add_motif(nway);
        mt2->add_motif(m2);
        mt2->add_motif(m3, 1);
        
        auto mset2 = motif_data_structure::MotifStateEnsembleTree(mt2);
        
        CHECK(mset2.size() == 4);
    }

    // graph tests!
    SUBCASE("test adding ensembles to graph") {

        auto mseg = motif_data_structure::MotifStateEnsembleGraph();
        auto mse = resources::Manager::instance().motif_state_ensemble("GG_LL_CC_RR");
        mseg.add_ensemble(*mse);
        mseg.add_ensemble(*mse);
        CHECK(mseg.size() == 2);

        auto mseg2 = motif_data_structure::MotifStateEnsembleOPGraph();
        mseg2.add_ensemble(mse);
        mseg2.add_ensemble(mse);
        CHECK(mseg2.size() == 2);
        CHECK(mseg2.get_ensemble(0) == mseg2.get_ensemble(1));

    }

    SUBCASE("setup from motif tree with graph") {
        /*auto mg = std::make_shared<motif_data_structure::MotifGraph>();
        auto m1 = resources::Manager::instance().motif("HELIX.IDEAL");
        auto m2 = resources::Manager::instance().motif("HELIX.IDEAL");
        mg->add_motif(m1);
        mg->add_motif(m2);*/
    }
}