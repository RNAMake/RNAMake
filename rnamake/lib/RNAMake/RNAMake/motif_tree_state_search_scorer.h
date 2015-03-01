//
//  motif_tree_state_search_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_search_scorer__
#define __RNAMake__motif_tree_state_search_scorer__

#include <stdio.h>
#include "motif_tree_state_search_node.h"
#include "basepair_state.h"

class MotifTreeStateSearchScorer {
public:
    MotifTreeStateSearchScorer() {}
    
    ~MotifTreeStateSearchScorer() {}
    
public:
    
    void
    setup(BasepairStateOP const &);
    
    virtual
    inline
    float
    score(MotifTreeStateSearchNodeOP const &) { return -1000; }
    
    float
    accept_score(MotifTreeStateSearchNodeOP const &);
    
protected:
    BasepairStateOP target_, target_flip_;
    
};

class MTSS_GreedyBestFirstSearch : public MotifTreeStateSearchScorer {
public:
    MTSS_GreedyBestFirstSearch() {}

public:

    inline
    float
    score(MotifTreeStateSearchNodeOP const & node) {
        return new_score_function(node->active_states()[0], target_, target_flip_);
    }
};

class MTSS_Astar: public MotifTreeStateSearchScorer {
public:
    MTSS_Astar() {
        g_ = 0; h_ = 0;
        ss_score_weight_ = 0.25;
        level_weight_ = 2.0;
    }
    
    inline
    float
    score(MotifTreeStateSearchNodeOP const & node) {
        g_ = node->ss_score()*ss_score_weight_;
        if(node->level() > 2) {
            g_ += node->level()*level_weight_;
        }
        h_ = new_score_function(node->active_states()[0], target_, target_flip_);
        return g_ + h_;
    }

private:
    float g_, h_;
    float ss_score_weight_, level_weight_;
};

typedef std::shared_ptr<MotifTreeStateSearchScorer> MotifTreeStateSearchScorerOP;


#endif /* defined(__RNAMake__motif_tree_state_search_scorer__) */
