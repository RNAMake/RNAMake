//
//  thermo_fluc_sampler.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluc_sampler.h"

void
ThermoFlucSampler::setup(
    MotifStateEnsembleTreeOP const & mset) {
    mset_ = mset;
    mst_ = mset_->to_mst();
    //set MonteCarlo temperature to kBT, Boltzmann constant in pN.A/K
    mc_ = MonteCarlo(temperature_*1.3806488e-1);
    states_ = Ints(mset_->size());
    for(int i = 0; i < mset_->size(); i++) { states_[i] = 0; }
    
}

void
ThermoFlucSampler::sample(
    int steps) {
    for(int i = 0; i < steps; i++) { next(); }
}

int
ThermoFlucSampler::next() {
    
    node_num_ = 1 + rng_.randrange((int)mst_->size()-2);
    mset_node_ = mset_->get_node(node_num_);
    mst_node_  = mst_->get_node(node_num_);
    pos_ = states_[node_num_];

    energy_ = mset_node_->data()->get_member(pos_)->energy;
    new_mem_ = mset_node_->data()->get_random_member();
    
    accept_ = mc_.accept(energy_, new_mem_->energy);
    
    if(accept_) {
        update(node_num_, new_mem_);
        return 1;
    }
    else {
        return 0;
    }
}

void
ThermoFlucSampler::update(
    int node_num,
    MotifStateEnsembleMemberOP const & new_mem) {
    
    last_state_pos_ = states_[node_num];
    states_[node_num_] = mset_node_->data()->member_index(new_mem);
    last_state_ = mst_node_->data()->ref_state;
    last_num_ = node_num;
    mst_->replace_state(node_num, new_mem->motif_state);
}

void
ThermoFlucSampler::to_pdb(
    String fname) {
    mst_->to_motif_tree()->to_pdb(fname);
}