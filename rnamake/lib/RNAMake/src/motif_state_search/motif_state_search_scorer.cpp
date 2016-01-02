//
//  motif_state_search_scorer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_state_search/motif_state_search_scorer.h"

void
MotifStateSearchScorer::set_target(
    BasepairStateOP const & target)
{
    target_ = target;
    target_flip_ = BasepairStateOP ( new BasepairState ( target_->copy() ));
    target_flip_->flip();
    
}

float
MotifStateSearchScorer::accept_score(
    MotifStateSearchNodeOP const & node) {
    
    best_score_ = 1000;
    int i = -1;
    for(auto const & state : node->cur_state()->end_states() ) {
        i++;
        if(i == 0) { continue; }
        
        //score_ = state->d().distance(target_->d());
        
        score_ = (state->sugars()[0].distance(target_->sugars()[1]) +
                  state->sugars()[1].distance(target_->sugars()[0]))*0.50;
        
        r_diff_       = state->r().difference(target_->r());
        r_diff_flip_  = state->r().difference(target_flip_->r());
    
        if( r_diff_ > r_diff_flip_) { score_ += 2*r_diff_flip_; }
        else                        { score_ += 2*r_diff_;      }
        
        if(score_ < best_score_) {
            best_score_ = score_;
        }
    
    }
    
    return best_score_;
}

