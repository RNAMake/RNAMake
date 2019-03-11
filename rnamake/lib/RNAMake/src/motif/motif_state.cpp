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

String
MotifState::to_str() {
    String s = name_ + "|" + std::to_string(score_) + "|" + std::to_string(size_) + "|";
    s += std::to_string(block_end_add_) + "|" + vectors_to_str(beads_) + "|";
    s += base::join_by_delimiter(end_names_, ",") + "|" + base::join_by_delimiter(end_ids_, ",") + "|";
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
        if(n == name) { return end_states_[i]; }
    }
    
    throw MotifStateException(
        "cannot find end_name: " + name + " in MotifState::get_end_state");
    
}

int
MotifState::get_end_index(
    BasepairStateOP const & end_state) {
    
    if(std::find(end_states_.begin(), end_states_.end(), end_state) == end_states_.end()) {
        throw MotifStateException(
            "cannot find end_state in end_states_ in MotifState::end_index");
    }
    int pos = (int)(std::find(end_states_.begin(), end_states_.end(), end_state) - end_states_.begin());
    return pos;
}

int
MotifState::get_end_index(
    String const & name) {
    
    int i = -1;
    for(auto const & n : end_names_) {
        i++;
        if(n == name) { return i; }
    }
    
    throw MotifStateException(
        "cannot find end_name: " + name + " in MotifState::get_end_index");
}






























