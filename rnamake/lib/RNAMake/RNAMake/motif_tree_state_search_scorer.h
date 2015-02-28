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
    MotifTreeStateSearchScorer(
        BasepairStateOP const &);
    
    ~MotifTreeStateSearchScorer() {}
    
public:
    
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
    MTSS_GreedyBestFirstSearch(BasepairStateOP const & target) : MotifTreeStateSearchScorer(target) {}

public:

    inline
    float
    score(MotifTreeStateSearchNodeOP const & node) {
        return new_score_function(node->active_states()[0], target_, target_flip_);
    }
};

typedef std::shared_ptr<MotifTreeStateSearchScorer> MotifTreeStateSearchScorerOP;


#endif /* defined(__RNAMake__motif_tree_state_search_scorer__) */
