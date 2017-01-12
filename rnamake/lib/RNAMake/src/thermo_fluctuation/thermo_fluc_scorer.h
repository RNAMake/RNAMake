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
    FrameScorer() : ThermoFlucScorer(),
        frame_score_( 0 ),
        r_diff_( 0 ),
        r_diff_flip_( 0 )
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
        state_2->flip();
        
        if(r_diff_ > r_diff_flip_) { frame_score_ += r_diff_flip_; }
        else                       { frame_score_ += r_diff_;      }
        
        return frame_score_;
        
    }
    
private:
    float frame_score_ = 0, r_diff_ = 0, r_diff_flip_ = 0;
};


class FrameScorerDevel : public ThermoFlucScorer {
public:
    FrameScorerDevel() : ThermoFlucScorer(),
        frame_score_( 0 ),
        r_diff_( 0 ),
        r_diff_flip_( 0 ),
        weight_d_( 1 ),
        weight_r_( 1 )
    {}
    
    ~FrameScorerDevel() {}
    
public:
    
    inline
    float
    score(
        BasepairStateOP & state_1,
        BasepairStateOP & state_2) {
        
        frame_score_ = state_1->d().distance(state_2->d())*weight_d_;
        r_diff_ = state_1->r().difference(state_2->r());
        state_2->flip();
        r_diff_flip_ = state_1->r().difference(state_2->r());;
        state_2->flip();
        if(r_diff_ > r_diff_flip_) { frame_score_ += r_diff_flip_*weight_r_; }
        else                       { frame_score_ += r_diff_*weight_r_;      }
        
        return frame_score_;
        
    }
    
public:
    
    inline
    void
    weight_d(
             float weight_d) {  weight_d_ = weight_d; }
    
    inline
    void
    weight_r(
             float weight_r) {  weight_r_ = weight_r; }
    

public:
    
    inline
    float
    weight_d() { return weight_d_; }
    
    inline
    float
    weight_r() { return weight_r_; }
    
private:
    float frame_score_ = 0, r_diff_ = 0, r_diff_flip_ = 0;
    float weight_d_ = 0, weight_r_ = 0;
};



typedef std::shared_ptr<ThermoFlucScorer> ThermoFlucScorerOP;

#endif /* defined(__RNAMake__thermo_fluc_scorer__) */
