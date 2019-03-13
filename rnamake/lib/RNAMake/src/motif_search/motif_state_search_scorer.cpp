//
//  motif_state_search_scorer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_search/motif_state_search_scorer.h"

namespace motif_search {

void
MotifStateSearchScorer::set_target(
        structure::BasepairStateOP const & target,
        bool target_an_aligned_end) {
    target_ = target;
    target_an_aligned_end_ = target_an_aligned_end;
    target_flip_ = structure::BasepairStateOP(new structure::BasepairState(target_->copy()));
    target_flip_->flip();

}

float
MotifStateSearchScorer::accept_score(
        MotifStateSearchNodeOP const & node) {
    best_score_ = 1000;
    int i = -1;
    for (auto const & state : node->cur_state()->end_states()) {
        i++;
        if (i == 0) { continue; }

        score_ = state->d().distance(target_->d());

        if (target_an_aligned_end_) {
            r_diff_ = state->r().difference(target_->r());
        } else {
            r_diff_ = state->r().difference(target_flip_->r());
        }
        score_ += 2 * r_diff_;

        if (score_ < best_score_) {
            best_score_ = score_;
        }
    }
    return best_score_;
}

}