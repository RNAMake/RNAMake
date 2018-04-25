
//headers for testing
#include "../common.hpp"
#include "../tools/motif_tree_builder.hpp"

//RNAMake Headers
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "motif_data_structures/motif_state_ensemble_graph.h"



TEST_CASE( "Test Assembling MotifEnsembleStates together", "[MotifStateEnsembleTree]" ) {
    
    SECTION("test adding ensembles to tree") {
    
        auto mset = MotifStateEnsembleTree();
        auto mse = RM::instance().motif_state_ensemble("GG_LL_CC_RR");
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        
        REQUIRE(mset.size() == 2);
    }
    
    SECTION("test getting a motif state tree from the ensemble tree") {
        
        auto mset = MotifStateEnsembleTree();
        auto mse = RM::instance().motif_state_ensemble("GG_LL_CC_RR");
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);
        mset.add_ensemble(mse);

        auto mst = mset.to_mst();
        REQUIRE(mst->size() == mset.size());
        
    }
    
    SECTION("setup from motif tree") {
        auto mt = std::make_shared<MotifTree>();
        auto m1 = RM::instance().motif("HELIX.IDEAL");
        auto m2 = RM::instance().motif("HELIX.IDEAL");
        mt->add_motif(m1);
        mt->add_motif(m2);

        auto ms1 = RM::instance().motif_state("HELIX.IDEAL");
        
        auto mset = MotifStateEnsembleTree(mt);
        
        REQUIRE(mset.size() == 2);
        
        auto mst = mset.to_mst();
        //REQUIRE(mst->get_node(0)->data()->name() == "CG=CG.0");
        
        auto mt2 = std::make_shared<MotifTree>();
        auto m3 = RM::instance().motif("HELIX.IDEAL");
        auto nway = RM::instance().motif("NWAY.1GID.0");
        
        mt2->add_motif(m1);
        mt2->add_motif(nway);
        mt2->add_motif(m2);
        mt2->add_motif(m3, 1);
        
        auto mset2 = MotifStateEnsembleTree(mt2);
        
        REQUIRE(mset2.size() == 4);
    }

    // graph tests!
    SECTION("test adding ensembles to graph") {

        auto mseg = MotifStateEnsembleGraph();
        auto mse = RM::instance().motif_state_ensemble("GG_LL_CC_RR");
        mseg.add_ensemble(mse);
        mseg.add_ensemble(mse);

        REQUIRE(mseg.size() == 2);
    }

    SECTION("setup from motif tree with graph") {
        auto mg = std::make_shared<MotifGraph>();
        auto m1 = RM::instance().motif("HELIX.IDEAL");
        auto m2 = RM::instance().motif("HELIX.IDEAL");
        mg->add_motif(m1);
        mg->add_motif(m2);



    }
}