//
//  motif_state_search_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_unittest__
#define __RNAMake__motif_state_search_unittest__

#include <stdio.h>

#include "unittest.h"
#include "motif_state_search/motif_state_search.h"

class aligner {
public:
    aligner():
    bp_state_(BasepairState()),
    bp_state_final_(BasepairState())
    {}
    
    ~aligner() {}
    
public:
    
    inline
    void
    get_aligned_motif_state(
        BasepairStateOP const & ref_bp_state,
        MotifStateOP & cur_state,
        MotifStateOP const & org_state) {
        
        ref_bp_state->get_transforming_r_and_t(*org_state->end_states()[0], bp_state_);
        for(int i = 0; i < org_state->end_states().size(); i++) {
            org_state->end_states()[i]->get_transformed_state(bp_state_, bp_state_final_);
            cur_state->update_end_state(i, bp_state_final_);
        }
        
        t_beads_ =  Points(org_state->beads().size());
        dot_vectors(bp_state_.r_T(), org_state->beads(), t_beads_);
        for(int i = 0; i < t_beads_.size();  i++) { t_beads_[i] += bp_state_.d(); }
        cur_state->beads(t_beads_);
    }
    
private:
    BasepairState bp_state_, bp_state_final_;
    Points t_beads_;
};


class MotifStateSearchUnittest : public Unittest {
public:
    MotifStateSearchUnittest() {}
    
    ~MotifStateSearchUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_search();
    
    int
    test_aligner();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_state_search_unittest__) */
