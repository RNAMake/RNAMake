//
//  motif_state_aligner.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_aligner__
#define __RNAMake__motif_state_aligner__

#include <stdio.h>
#include "motif/motif_state.h"
#include "structure/basepair_state.h"

namespace motif {

class MotifStateAligner {
public:
    MotifStateAligner() {}

    ~MotifStateAligner() {}

public:
    inline
    void
    get_aligned_motif_state(
            structure::BasepairStateOP const & ref_bp_state,
            MotifStateOP & cur_state,
            MotifStateOP const & org_state) {

        ref_bp_state->get_transforming_r_and_t(*org_state->end_states()[0], bp_state_);
        for (int i = 0; i < org_state->end_states().size(); i++) {
            org_state->end_states()[i]->get_transformed_state(bp_state_, bp_state_final_);
            cur_state->update_end_state(i, bp_state_final_);
        }

        t_beads_ = math::Points(org_state->beads().size());
        math::dot_vectors(bp_state_.r_T(), org_state->beads(), t_beads_);
        for (int i = 0; i < t_beads_.size(); i++) { t_beads_[i] += bp_state_.d(); }
        cur_state->beads(t_beads_);
    }

private:
    structure::BasepairState bp_state_, bp_state_final_;
    math::Points t_beads_;

};

}

#endif /* defined(__RNAMake__motif_state_aligner__) */
