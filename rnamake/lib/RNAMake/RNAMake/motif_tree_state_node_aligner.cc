//
//  motif_tree_state_node_aligner.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_node_aligner.h"

void
MotifTreeStateNodeAligner::transform_state(
    BasepairStateOP const & parent_end,
    MotifTreeStateNodeOP const &parent,
    MotifTreeStateNodeOP const & child) {
    
    parent_end->get_transforming_r_and_t(ref_bp_state_, r_state_);
    for (auto const & s : child->states()) {
        if(s == NULL) { continue; }
        s->get_transformed_state(r_state_, t_state_);
        s->set(t_state_);
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
