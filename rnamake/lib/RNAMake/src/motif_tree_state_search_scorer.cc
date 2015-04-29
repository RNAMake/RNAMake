//
//  motif_tree_state_search_scorer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_search_scorer.h"

void
MotifTreeStateSearchScorer::setup(
    BasepairStateOP const & target)
{
    target_ = target;
    target_flip_ = BasepairStateOP ( new BasepairState ( target_->copy() ));
    target_flip_->flip();
    
}

float
MotifTreeStateSearchScorer::accept_score(
    MotifTreeStateSearchNodeOP const & node) {
    
    BasepairStateOP current = node->active_states()[0];
    float score = current->d().distance(target_->d());
    float r_diff       = current->r().difference(target_->r());
    float r_diff_flip  = current->r().difference(target_flip_->r());
    
    if( r_diff > r_diff_flip) {
        r_diff = r_diff_flip;
    }
    score += 2*r_diff;
    return score;
}






