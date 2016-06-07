//
//  motif_state.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdexcept>
#include <algorithm>

#include "math/xyz_vector.h"
#include "structure/basepair_state.h"
#include "motif/motif_state.h"

MotifState
MotifState::copy() {
    BasepairStateOPs end_states;
    for(auto const & bp_state : end_states_) {
        end_states.push_back(std::make_shared<BasepairState>(bp_state->copy()));
    }
    Points beads;
    for (auto const & b : beads_) {
        beads.push_back(Point(b));
    }
    
    return MotifState(name_, end_names_, end_ids_, end_states, beads, score_,
                    size_, block_end_add_);
}

String
MotifState::to_str() {
    String s = name_ + "|" + std::to_string(score_) + "|" + std::to_string(size_) + "|";
    s += std::to_string(block_end_add_) + "|" + vectors_to_str(beads_) + "|";
    s += join_by_delimiter(end_names_, ",") + "|" + join_by_delimiter(end_ids_, ",") + "|";
    for(auto const & state : end_states_) {
        s += state->to_str() + "|";
    }
    return s;
}


BasepairStateOP const &
MotifState::get_end_state(
    String const & name) {
    
    int i = -1;
    for(auto const & n : end_names_) {
        i++;
        if(i == block_end_add_) { continue; } 
        if(n == name) { return end_states_[i]; }
    }
    
    throw std::runtime_error("cannot find end_name: " + name + " in MotifState::get_end_state");
    
}

int
MotifState::end_index(
    BasepairStateOP const & end_state) {
    if(std::find(end_states_.begin(), end_states_.end(), end_state) == end_states_.end()) {
        throw std::runtime_error("cannot find end_state in end_states_ in MotifState::end_index");
    }
    int pos = (int)(std::find(end_states_.begin(), end_states_.end(), end_state) - end_states_.begin());
    return pos;
}































