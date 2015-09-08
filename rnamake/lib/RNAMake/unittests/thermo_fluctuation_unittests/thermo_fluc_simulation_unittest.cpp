//
//  thermo_fluc_simulation_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluc_simulation_unittest.h"
#include "thermo_fluctuation/thermo_fluc_simulation.h"


int
ThermoFlucSimulationUnittest::test_creation() {
    ThermoFlucSimulation tfs;
    tfs.option("steps", 10);
    if(tfs.option<int>("steps") != 10) { return 0; }
    return 1;
}

int
ThermoFlucSimulationUnittest::test_run() {
    return 1;
}

int
ThermoFlucSimulationUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    if (test_run() == 0)       {  std::cout << "test_run failed" << std::endl; }
    return 1;
}