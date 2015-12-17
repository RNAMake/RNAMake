//
//  motif_state_search_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_scorer__
#define __RNAMake__motif_state_search_scorer__

#include <stdio.h>

//RNAMAke Headers
#include "structure/basepair_state.h"
#include "motif_state_search/motif_state_search_node.h"

class MotifStateSearchScorer {
public:
    MotifStateSearchScorer() {}
    
    ~MotifStateSearchScorer() {}
    
public:
    
    void
    set_target(BasepairStateOP const &);
    
    virtual
    inline
    float
    score(
        MotifStateSearchNodeOP const & node ) {
        best_score_ = 1000;
        int i = -1;
        for(auto const & state : node->cur_state()->end_states() ) {
            i++;
            if(i == 0) { continue; }
            
            score_ = new_score_function_new(state, target_, target_flip_);
            
            if(score_ < best_score_) {
                best_score_ = score_;
            }
            
        }
        return best_score_;
    
    }
    
    float
    accept_score(MotifStateSearchNodeOP const &);
    
public:
    virtual
    inline
    void
    level_weight(float const & nlevel_weight) {}
    
    virtual
    inline
    void
    ss_score_weight(float const & nss_score_weight) {}
    
protected:
    BasepairStateOP target_, target_flip_;
    float best_score_, score_, r_diff_, r_diff_flip_;
    
};

class MSS_GreedyBestFirstSearch : public MotifStateSearchScorer {
public:
    MSS_GreedyBestFirstSearch() {}
    
public:
    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        return MotifStateSearchScorer::score(node);
    }
    
public:
    
};

class MSS_Astar: public MotifStateSearchScorer {
public:
    MSS_Astar() {
        g_ = 0; h_ = 0;
        ss_score_weight_ = 0.25;
        level_weight_ = 2.0;
    }
    
    inline
    float
    score(MotifStateSearchNodeOP const & node) {
        g_ = node->ss_score()*ss_score_weight_;
        if(node->level() > 2) {
            g_ += node->level()*level_weight_;
        }
        h_ = MotifStateSearchScorer::score(node);
        return g_ + h_;
    }
    
public:
    
    inline
    void
    level_weight(float const & nlevel_weight) { level_weight_ = nlevel_weight; }
    
    inline
    void
    ss_score_weight(float const & nss_score_weight) { ss_score_weight_ = nss_score_weight; }
    
private:
    float g_, h_;
    float ss_score_weight_, level_weight_;
};

typedef std::shared_ptr<MotifStateSearchScorer> MotifStateSearchScorerOP;

#endif /* defined(__RNAMake__motif_state_search_scorer__) */
