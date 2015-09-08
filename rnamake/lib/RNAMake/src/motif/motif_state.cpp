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


MotifState
str_to_motif_state(
    String const & s) {
    
    auto spl = split_str_by_delimiter(s, "|");
    String name = spl[0];
    float score = std::stof(spl[1]);
    int size = std::stoi(spl[2]);
    int block_end_add = std::stoi(spl[3]);
    Points beads = vectors_from_str(spl[4]);
    auto end_names = split_str_by_delimiter(spl[5], ",");
    auto end_ids = split_str_by_delimiter(spl[6], ",");
    BasepairStateOPs end_states;
    for(int i = 7; i < spl.size(); i++) {
        auto end_state = std::make_shared<BasepairState>(str_to_basepairstate(spl[i]));
        end_states.push_back(end_state);
    }
    return MotifState(name, end_names, end_ids, end_states, beads, score, size, block_end_add);
    
}

































