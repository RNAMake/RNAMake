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
    
    void
    transform_state(
        BasepairStateOP,
        MotifTreeStateNodeOP ,
        MotifTreeStateNodeOP);

    
    void
    transform_beads(
        MotifTreeStateNodeOP const &);
    
    
private:
    BasepairState r_state_, t_state_, ref_bp_state_;
    
};


#endif /* defined(__RNAMake__motif_tree_state_node_aligner__) */
