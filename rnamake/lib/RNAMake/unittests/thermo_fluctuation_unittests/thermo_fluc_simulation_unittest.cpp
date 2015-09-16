//
//  thermo_fluc_simulation_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluc_simulation_unittest.h"
#include "thermo_fluctuation/thermo_fluc_simulation.h"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"


int
ThermoFlucSimulationUnittest::test_creation() {
    ThermoFlucSimulation tfs;
    tfs.option("steps", 10);
    if(tfs.option<int>("steps") != 10) { return 0; }
    return 1;
}

int
ThermoFlucSimulationUnittest::test_run() {
    ThermoFlucSimulation tfs;
    auto mset = std::make_shared<MotifStateEnsembleTree>();
    auto mse = ResourceManager::getInstance().get_motif_state_ensemble("GG_LL_CC_RR");
    for(int i = 0; i < 10; i++) { mset->add_ensemble(mse); }
    tfs.setup(mset, 0, 9, 0, 1);
    tfs.option("cutoff", 30.0f);
    int count = tfs.run();    
    return 1;
}

int
ThermoFlucSimulationUnittest::run() {
    //if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    //if (test_run() == 0)       {  std::cout << "test_run failed" << std::endl; }
    return 1;
}