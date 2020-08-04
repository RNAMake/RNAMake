//
//  thermo_fluc_sampler.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "thermo_fluctuation/thermo_fluc_sampler.h"

namespace thermo_fluctuation {

void
ThermoFlucSampler::setup(
        motif_data_structure::MotifStateEnsembleTreeOP const & mset) {
    mset_ = mset;
    mst_ = mset_->to_mst();
    //set MonteCarlo temperature to kBT, Boltzmann constant in pN.A/K
    mc_ = util::MonteCarlo(temperature_);
    states_ = Ints(mset_->size());
    for (int i = 0; i < mset_->size(); i++) { states_[i] = 0; }

    if(randomized_start_) {
        for(int i = 0; i < mset_->size(); i++) {
            node_num_ = i;
            mset_node_ = mset_->get_node(i);
            mst_node_ = mst_->get_node(i);
            mem_pos_ = rng_.randrange(mset_node_->data()->size());
            if (mem_pos_ == mset_node_->data()->size()) { mem_pos_--; }
            new_mem_ = mset_node_->data()->get_member(mem_pos_);
            update(i, new_mem_);
        }
    }

}

void
ThermoFlucSampler::sample(
        int steps) {
    for (int i = 0; i < steps; i++) { next(); }
}

int
ThermoFlucSampler::next() {

    node_num_ = rng_.randrange((int) mst_->size());
    if (node_num_ == mst_->size()) { node_num_--; }

    mset_node_ = mset_->get_node(node_num_);
    mst_node_ = mst_->get_node(node_num_);
    pos_ = states_[node_num_];
    mem_pos_ = rng_.randrange(mset_node_->data()->size());
    if (mem_pos_ == mset_node_->data()->size()) { mem_pos_--; }

    energy_ = mset_node_->data()->get_member(pos_)->energy;
    new_mem_ = mset_node_->data()->get_member(mem_pos_);

    accept_ = mc_.accept(energy_, new_mem_->energy);

    if (accept_) {
        update(node_num_, new_mem_);
        return 1;
    } else {
        return 0;
    }
}

void
ThermoFlucSampler::update(
        int node_num,
        motif::MotifStateEnsembleMemberOP const & new_mem) {

    last_state_pos_ = states_[node_num];
    states_[node_num_] = mset_node_->data()->member_index(new_mem);
    last_state_ = mst_node_->data()->ref_state;
    last_num_ = node_num;
    mst_->replace_state(node_num, new_mem->motif_state);
}

void
ThermoFlucSampler::to_pdb(
        String fname,
        int renumber) {
    mst_->to_motif_tree()->to_pdb(fname, renumber);
}

}