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

Ints const
MotifTreeStateNode::available_ends() {
    Ints indices;
    int i = -1;
    for( auto const & state : states_) {
        i++;
        if(state == NULL) { continue; }
        if(children_[i] == NULL) {
            indices.push_back(i);
        }
    }
    return indices;
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
MotifTreeStateNode::replace_mts(MotifTreeStateOP const & nmts) {
    mts_ = nmts;
    int i = -1;
    for (auto const & s : mts_->end_states()) {
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
    int const & i) {
    children_[i] = node;
}

MotifTreeStateNodeOPs
MotifTreeStateNode::children() {
    MotifTreeStateNodeOPs children;
    for(auto const & c : children_) {
        if(c != NULL) { children.push_back(c); }
    }
    return children;
}




