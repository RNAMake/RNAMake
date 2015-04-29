//
//  motif_tree_state_search_node.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/23/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_node.h"
#include "motif_tree_state_search_node.h"


MotifTreeStateSearchNode::MotifTreeStateSearchNode(
    MotifTreeStateOP const & mts,
    MotifTreeStateSearchNodeOP const & parent,
    int lib_type):
    mts_ ( mts ),
    parent_ ( parent ),
    lib_type_ (lib_type )
{
        
    if(parent == NULL) { level_ = 0; }
    else               { level_ = parent->level() + 1; }
    
    score_ = 1000; size_ = 0; ss_score_ = 0.0f;
    beads_ = Points();
    node_counts_ = Ints();
    states_ = BasepairStateOPs(mts->end_states().size());
    active_ = Ints(mts->end_states().size());
    int i = 0;
    for (auto const & s : mts_->end_states()) {
        if( s == NULL) {
            states_[i] = BasepairStateOP(new BasepairState());
            active_[i] = 0;
        }
        else           {
            states_[i] = BasepairStateOP(new BasepairState(s->copy()));
            active_[i] = 1;
        }
        i++;
    }
        
}


void
MotifTreeStateSearchNode::replace_mts(MotifTreeStateOP const & nmts) {
    mts_ = nmts;
    int i = -1;
    if(states_.size() < mts_->end_states().size()) {
        states_.push_back(BasepairStateOP(new BasepairState()));
        active_.push_back(0);
    }
    
    for (auto const & s : mts_->end_states()) {
        i++;
        if(s == NULL) { active_[i] = 0; }
        else          { active_[i] = 1; }
    }
}


BasepairStateOPs
MotifTreeStateSearchNode::active_states() const {
    BasepairStateOPs states;
    int i = -1;
    for (auto const & s : states_ ) {
        i++;
        if(active_[i] == 1) { states.push_back(s); }
    }
    return states;
}

