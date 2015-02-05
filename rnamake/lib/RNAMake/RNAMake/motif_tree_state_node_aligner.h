//
//  motif_tree_state_node_aligner.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_node_aligner__
#define __RNAMake__motif_tree_state_node_aligner__

#include <stdio.h>
#include "xyzMatrix.h"
#include "basepair_state.h"
#include "motif_tree_state_node.fwd.h"
#include "motif_tree_state_node.h"

class MotifTreeStateNodeAligner {
public:
    MotifTreeStateNodeAligner():
    r_state_ ( BasepairState() ),
    t_state_ ( BasepairState() ),
    ref_bp_state_ ( get_ref_bp_state() )
    {}
    
    ~MotifTreeStateNodeAligner() {}
    
public:
    
    inline
    void
    transform_state(
        BasepairStateOP const & parent_end,
        MotifTreeStateNodeOP const & parent,
        MotifTreeStateNodeOP const & child) {
        
        parent_end->get_transforming_r_and_t(ref_bp_state_, r_state_);
        for (const auto & s : child->states()) {
            if(s == NULL) { continue; }
            s->get_transformed_state(r_state_, t_state_);
            s->set(t_state_);
        }

    }
    
    inline
    void
    transform_beads(
        MotifTreeStateNodeOP const & child) {
        Points t_beads (child->mts().beads().size());
        dot_vectors(r_state_.r_T(), child->mts().beads(), t_beads);
        for(int i = 0; i < t_beads.size();  i++) { t_beads[i] += r_state_.d(); }
        child->beads(t_beads);
    }
    
    
private:
    BasepairState r_state_, t_state_, ref_bp_state_;
    
};


#endif /* defined(__RNAMake__motif_tree_state_node_aligner__) */
