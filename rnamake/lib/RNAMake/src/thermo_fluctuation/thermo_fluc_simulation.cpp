//
//  thermo_fluc_simulation.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluctuation/thermo_fluc_simulation.h"


void
ThermoFlucSimulation::setup_options() {
    options_ = Options();
    options_.add_option(Option("temperature", 298.15f));
    options_.add_option(Option("steps", 100000));
    options_.add_option(Option("record", 1));
    options_.add_option(Option("cutoff", 5.0f));
    update_var_options();
}

void
ThermoFlucSimulation::update_var_options() {
    temperature_   = options_.option<float>("temperature");
    steps_         = options_.option<int>("steps");
    record_        = options_.option<int>("record");
    cutoff_        = options_.option<float>("cutoff");
}

void
ThermoFlucSimulation::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    sampler_.setup(mset);
    sampler_.temperature(option<float>("temperature"));
    
    n1_ = sampler_.mst()->get_node(ni1);
    n2_ = sampler_.mst()->get_node(ni2);
    ei1_ = ei1;
    ei2_ = ei2;
}

int
ThermoFlucSimulation::run() {
    
    int steps = 0;
    int r = 0;
    while (steps < steps_) {
        r = sampler_.next();
        
        std::cout << r << std::endl;
        
        steps++;
    }
    
    return 0;
}