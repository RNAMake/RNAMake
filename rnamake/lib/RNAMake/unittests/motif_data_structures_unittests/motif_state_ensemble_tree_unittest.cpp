//
//  motif_state_ensemble_tree_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "build_motif_tree.h"
#include "motif_state_ensemble_tree_unittest.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "resources/resource_manager.h"
#include "resources/motif_state_ensemble_sqlite_library.h"


int
MotifStateEnsembleTreeUnittest::test_creation() {
    MotifStateEnsembleTree mse_tree;
    return 1;
}

int
MotifStateEnsembleTreeUnittest::test_add_ensemble() {
    MotifStateEnsembleTree mse_tree;
    auto mse = ResourceManager::getInstance().get_motif_state_ensemble("GG_LL_CC_RR");
    mse_tree.add_ensemble(mse);
    mse_tree.add_ensemble(mse);
    if(mse_tree.size() != 2) { return 0; }
    return 1;
}

int
MotifStateEnsembleTreeUnittest::test_to_mst() {
    MotifStateEnsembleTree mset;
    auto mse = ResourceManager::getInstance().get_motif_state_ensemble("GG_LL_CC_RR");
    for(int i = 0; i < 10; i++) { mset.add_ensemble(mse); }
    
    auto mst = mset.to_mst();
    mst->to_motif_tree()->write_pdbs();
    return 1;
    
}

int
MotifStateEnsembleTreeUnittest::test_from_mt() {
    BuildMotifTree builder(Strings{"bp_steps", "twoway"});
    auto mt = builder.build(10);
    MotifStateEnsembleTree mset;
    mset.setup_from_mt(mt);
    auto mst = mset.to_mst();
    
    return 1;
}

int
MotifStateEnsembleTreeUnittest::test_enumerator() {
    MotifStateEnsembleSqliteLibrary lib("bp_steps");
    lib.load_all();
    for(auto const & mes : lib) {
    
        auto mtst = std::make_shared<MotifStateEnsembleTree>();
        mtst->add_ensemble(mes);
        MotifStateEnsembleTreeEnumerator enumerator(mtst);
        enumerator.record(mes->id());
    }
    
    return 1;
}


int
MotifStateEnsembleTreeUnittest::run() {
    //if (test_creation() == 0)      { std::cout << "test_creation failed" << std::endl;  }
    //if (test_add_ensemble() == 0)  { std::cout << "test_add_ensemble failed" << std::endl;  }
    //if (test_to_mst() == 0)        { std::cout << "test_to_mst failed" << std::endl;  }
    //if (test_from_mt() == 0)       { std::cout << "test_from_mt failed" << std::endl;  }
    if (test_enumerator() == 0)      { std::cout << "test_enumerator failed" << std::endl;  }

    return 1;
}