//
//  motif_ensemble.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_ensemble.h"

MotifState const &
MotifEnsemble::get_state(String const & name) const {
    for( auto const & me : motif_states_) {
        if(me.mts->name().compare(name) == 0) {
            return me;
        }
    }
    throw "cannot find mts";
}

MotifState const &
MotifEnsemble::get_random_state() const {
    int node_num = rand() % motif_states_.size();
    return motif_states_[node_num];
}