//
//  thermo_fluc_sampler_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluc_sampler_unittest.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"

int
ThermoFlucSamplerUnittest::test_creation() {
    ThermoFlucSampler tfs;
    return 1;
}

int
ThermoFlucSamplerUnittest::test_sample() {
    ThermoFlucSampler tfs;
    auto mset = std::make_shared<MotifStateEnsembleTree>();
    auto mse = ResourceManager::getInstance().get_motif_state_ensemble("GG_LL_CC_RR");
    for(int i = 0; i < 10; i++) { mset->add_ensemble(mse); }
    tfs.setup(mset);
    tfs.to_pdb("start.pdb");
    tfs.sample(10000);
    tfs.to_pdb("end.pdb");

    
    
    return 1;
}

int
ThermoFlucSamplerUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    if (test_sample() == 0)    {  std::cout << "test_sample failed" << std::endl; }
    return 0;
}