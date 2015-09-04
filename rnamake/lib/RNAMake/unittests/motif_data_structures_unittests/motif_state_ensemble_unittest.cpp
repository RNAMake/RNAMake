//
//  motif_state_ensemble_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_ensemble_unittest.h"
#include "resources/resource_manager.h"

int
MotifStateEnsembleUnittest::test_creation() {
    auto mse = ResourceManager::getInstance().get_motif_state_ensemble("GG_LL_CC_RR");
    auto mem = mse->get_random_member();
    std::cout << mem->energy << std::endl;
    return 1;
}

int
MotifStateEnsembleUnittest::run() {
    if (test_creation() == 0)      { std::cout << "test_creation failed" << std::endl;  }
    
    return 1;
}