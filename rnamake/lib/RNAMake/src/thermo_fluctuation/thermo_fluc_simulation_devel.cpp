//
//  thermo_fluc_simulation_devel.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//
//

#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"

void
ThermoFlucSimulationDevel::setup_options() {
    options_ = Options();
    options_.add_option(Option("temperature", 298.15f));
    options_.add_option(Option("steps", 100000));
    options_.add_option(Option("record", 1));
    options_.add_option(Option("cutoff", 5.0f));
    update_var_options();
}

void
ThermoFlucSimulationDevel::update_var_options() {
    temperature_   = options_.option<float>("temperature");
    steps_         = options_.option<int>("steps");
    record_        = options_.option<int>("record");
    cutoff_        = options_.option<float>("cutoff");
}

void
ThermoFlucSimulationDevel::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    sampler_.setup(mset);
    sampler_.temperature(option<float>("temperature"));
    
    ni1_ = ni1;
    ni2_ = ni2;
    ei1_ = ei1;
    ei2_ = ei2;
}

int
ThermoFlucSimulationDevel::run() {
    
    int steps = 0;
    int r = 0;
    int count = 0;
    int clash = 0;
    Ints check_nodes = { 22, 21 };
    Ints check_nodes_2 = { 1 };
    
    //fixes if length of tecto changes, need to come up with better system!
    check_nodes[0] = sampler_.mst()->last_node()->index();
    check_nodes[1] = sampler_.mst()->last_node()->index()-1;
    
    while (steps < steps_) {
        r = sampler_.next();
        //if(r == 0) { continue; }
        
        clash = 0;
        for(auto const & i : check_nodes) {
            for(auto const & j : check_nodes_2) {
                for(auto const & b2 : sampler_.mst()->get_node(i)->data()->cur_state->beads()) {
                    for(auto const & b1 : sampler_.mst()->get_node(j)->data()->cur_state->beads()) {
                        if(b1.distance(b2) < 2.2) { clash = 1; }
                    }
                    
                    if(clash) { break; }
                }
                if(clash) { break; }
            }
            if(clash) { break; }
        }
        
        if(clash) {
            //steps++;
            continue;
        }
        
        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->cur_state->end_states()[ei1_];
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->cur_state->end_states()[ei2_];
        
        score_ = scorer_->score(end_state_1_, end_state_2_);
        if(score_ < cutoff_) {
            /*try {
             sampler_.mst()->to_motif_tree()->to_pdb("test." + std::to_string(count) + ".pdb");
             }
             catch(...) { }*/
            count += 1;
        }
        
        steps++;
    }
    
    return count;
}