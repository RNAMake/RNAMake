//
//  motif_tree_state_node.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_node.h"

MotifTreeStateNode
MotifTreeStateNode::copy() {
    MotifTreeStateNode cmtn ( mts_, index_, parent_, lib_type_);
    cmtn.beads_ = beads_;
    cmtn.score_ = score_;
    cmtn.states_ = BasepairStateOPs(states_.size());
    int i = 0;
    for (auto const & s : states_) {
        if( s == NULL) { states_[i] = NULL; }
        else           {
            cmtn.states_[i] = BasepairStateOP(new BasepairState(s->copy()));
        }
        i++;
    }
    cmtn.children_ = MotifTreeStateNodeOPs(states_.size());
    return cmtn;
}

BasepairStateOPs const
MotifTreeStateNode::available_ends() {
    BasepairStateOPs states;
    int i = -1;
    for( auto const & state : states_) {
        i++;
        if(state == NULL) { continue; }
        if(children_[i] == NULL) {
            states.push_back(state);
        }
    }
    return states;
}

int
MotifTreeStateNode::parent_end_index() {
    if ( parent_ == NULL) { return -1; }
    int i = -1;
    for (auto const & n : parent_->children_) {
        i++;
        if(n.get() == this) {
            return i;
        }
    }
    throw "cannot find parent_end_index";
}

BasepairStateOP const &
MotifTreeStateNode::parent_end() {
    if ( parent_ == NULL) { return NULL; }
    int i = -1;
    for( auto const & n : parent_->children_) {
        i++;
        if(n.get() == this) {
            return parent_->states_[i];
        }
    }
    throw "cannot find parent_end";
}

void
MotifTreeStateNode::replace_mts(MotifTreeState const & nmts) {
    mts_ = nmts;
    int i = -1;
    for (auto const & s : mts_.end_states()) {
        i++;
        if(s == NULL) { states_[i] = NULL; }
        else {
            states_[i] = BasepairStateOP(new BasepairState(s->copy()));
        }
    }
}

void
MotifTreeStateNode::add_child(
    MotifTreeStateNodeOP const & node,
    BasepairStateOP const & end) {
    int i = (int)(std::find(states_.begin(), states_.end(), end) - states_.begin());
    children_[i] = node;
}





