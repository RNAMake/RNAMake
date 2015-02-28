//
//  motif_tree_state_node_aligner.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_node_aligner.h"
#include "motif_tree_state_search_node.h"

void
MotifTreeStateNodeAligner::transform_state(
    BasepairStateOP const & parent_end,
    MotifTreeStateNodeOP const &parent,
    MotifTreeStateNodeOP const & child) {
    
    parent_end->get_transforming_r_and_t(ref_bp_state_, r_state_);
    int i = -1;
    for (auto const & s : child->mts()->end_states()) {
        i++;
        if(s == NULL) { continue; }
        s->get_transformed_state(r_state_, t_state_);
        child->states()[i]->set(t_state_);
    }
    
}

void
MotifTreeStateNodeAligner::transform_beads(
    MotifTreeStateNodeOP const & child) {
    Points t_beads (child->mts()->beads().size());
    dot_vectors(r_state_.r_T(), child->mts()->beads(), t_beads);
    for(int i = 0; i < t_beads.size();  i++) { t_beads[i] += r_state_.d(); }
    child->beads(t_beads);
}


void
MotifTreeStateNodeAligner::transform_state(
    BasepairStateOP const & parent_end,
    MotifTreeStateSearchNodeOP const & parent,
    MotifTreeStateSearchNodeOP const & child) {
    
    parent_end->get_transforming_r_and_t(ref_bp_state_, r_state_);
    int i = -1;
    for (auto const & s : child->mts()->end_states()) {
        i++;
        if(s == NULL) { continue; }
        s->get_transformed_state(r_state_, t_state_);
        child->states()[i]->set(t_state_);
    }
    
}

void
MotifTreeStateNodeAligner::transform_beads(
    MotifTreeStateSearchNodeOP const & child) {
    Points t_beads (child->mts()->beads().size());
    dot_vectors(r_state_.r_T(), child->mts()->beads(), t_beads);
    for(int i = 0; i < t_beads.size();  i++) { t_beads[i] += r_state_.d(); }
    child->beads(t_beads);
}
