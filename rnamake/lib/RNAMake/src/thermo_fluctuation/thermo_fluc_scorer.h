//
//  thermo_fluc_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_scorer__
#define __RNAMake__thermo_fluc_scorer__

#include <stdio.h>

//RNAMake Headers
#include "structure/basepair_state.fwd.h"
#include "structure/basepair_state.h"


class ThermoFlucScorer {
public:
    ThermoFlucScorer() {}
    
    virtual
    ~ThermoFlucScorer() {}
    
    virtual
    inline
    float
    score(
        BasepairStateOP & state_1,
        BasepairStateOP & state_2) {
        return 0;
    }
    
};


class FrameScorer : public ThermoFlucScorer {
public:
    FrameScorer() : ThermoFlucScorer()
    {}
    
    ~FrameScorer() {}

public:
    
    inline
    float
    score(
        BasepairStateOP & state_1,
        BasepairStateOP & state_2) {
        
        frame_score_ = state_1->d().distance(state_2->d());
        r_diff_ = state_1->r().difference(state_2->r());
        state_2->flip();
        r_diff_flip_ = state_1->r().difference(state_2->r());;
        state_1->flip();
        
        if(r_diff_ > r_diff_flip_) { frame_score_ += r_diff_flip_; }
        else                       { frame_score_ += r_diff_;      }
        
        return frame_score_;
        
    }
    
private:
    float frame_score_, r_diff_, r_diff_flip_;
};


typedef std::shared_ptr<ThermoFlucScorer> ThermoFlucScorerOP;

#endif /* defined(__RNAMake__thermo_fluc_scorer__) */