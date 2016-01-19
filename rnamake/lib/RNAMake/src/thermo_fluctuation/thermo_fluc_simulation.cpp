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
    options_.add_option("temperature", 298.15f, OptionType::FLOAT);
    options_.add_option("steps", 100000, OptionType::INT);
    options_.add_option("record", false, OptionType::BOOL);
    options_.add_option("cutoff", 5, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
ThermoFlucSimulation::update_var_options() {
    temperature_   = options_.get_float("temperature");
    steps_         = options_.get_int("steps");
    cutoff_        = options_.get_float("cutoff");
    record_        = options_.get_bool("record");
}

void
ThermoFlucSimulation::setup(
    MotifStateEnsembleTreeOP const & mset,
    int ni1,
    int ni2,
    int ei1,
    int ei2) {
    
    sampler_.setup(mset);
    sampler_.temperature(temperature_);
    
    ni1_ = ni1;
    ni2_ = ni2;
    ei1_ = ei1;
    ei2_ = ei2;
}

int
ThermoFlucSimulation::run() {
    
    int steps = 0;
    int r = 0;
    int count = 0;
    int clash = 0;
    Ints check_nodes = { 22, 21 };
    Ints check_nodes_2 = { 1 };

    while (steps < steps_) {
        //if(r == 0) { continue; }
        
        end_state_1_ = sampler_.mst()->get_node(ni1_)->data()->cur_state->end_states()[ei1_];
        end_state_2_ = sampler_.mst()->get_node(ni2_)->data()->cur_state->end_states()[ei2_];
        
        int frame_score = end_state_1_->d().distance(end_state_2_->d());
        int r_diff = end_state_1_->r().difference(end_state_2_->r());
        end_state_2_->flip();
        int r_diff_flip = end_state_1_->r().difference(end_state_2_->r());;
        end_state_2_->flip();
        
        if(r_diff > r_diff_flip) {
            
        }
        

        std::cout << score_ << std::endl;
        
        exit(0);
        
        return score_;
       
        r = sampler_.next();
 
        
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