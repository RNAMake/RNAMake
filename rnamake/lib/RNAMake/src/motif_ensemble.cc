//
//  motif_ensemble.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


inline float GenXORShift32(void)
{
    static unsigned seed = 2463534242U;
    
    seed ^= (seed << 5);
    seed ^= (seed >> 13);
    seed ^= (seed << 6);
    
    return seed * (1.0f / 4294967295.0f);
}


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
    //int node_num = (int)((motif_states_.size()-1)*(rand() / (RAND_MAX + 1.0)));
    //int node_num = (int)((motif_states_.size())*dist_(mt_));
    int node_num = 1;
    
    return motif_states_[node_num];
}