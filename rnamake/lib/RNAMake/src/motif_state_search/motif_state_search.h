//
//  motif_state_search.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search__
#define __RNAMake__motif_state_search__

#include <stdio.h>

//RNAMake Headers
#include "base/base.h"
#include "base/option.h"
#include "motif/motif_state_aligner.h"
#include "motif_state_search/motif_state_selector.h"
#include "motif_state_search/motif_state_search_scorer.h"
#include "motif_state_search/motif_state_search_node.fwd.h"
#include "motif_state_search/motif_state_search_node.h"
#include "motif_state_search/motif_state_search_solution.h"


class MotifStateSearch : public Base {
public:
    MotifStateSearch():
    queue_(MotifStateSearchNodeQueue()),
    selector_(default_selector()),
    scorer_(std::make_shared<MotifStateSearchScorer>(MSS_Astar())),
    solutions_(MotifStateSearchSolutionOPs()),
    aligner_(MotifStateAligner())
    {
        setup_options();
    }
    
    ~MotifStateSearch() {}
    
public:
    
    MotifStateSearchSolutionOPs
    search(
        BasepairStateOP const &,
        BasepairStateOP const &);
    
    void
    setup(
        BasepairStateOP const &,
        BasepairStateOP const &);
    
    MotifStateSearchSolutionOP
    next();
    
    int
    finished();
    
public:
    inline
    void
    beads(Points const & beads) { beads_ = beads; }
    
protected:
    
    void
    setup_options();
    
    void
    update_var_options();
    
private:
    
    MotifStateSearchNodeOP
    _start_node(
        BasepairStateOP const &);
    
    MotifStateSearchSolutionOP
    _search();
    
    
private:
    MotifStateSearchNodeQueue queue_;
    MotifStateSelectorOP selector_;
    MotifStateSearchScorerOP scorer_;
    MotifStateSearchSolutionOPs solutions_;
    MotifStateSearchNodeOP test_node_;
    MotifStateAligner aligner_;
    MotifStateandTypes possible_children_;
    Points beads_;
    int no_more_solutions_;
    //options
    int sterics_, max_node_level_, min_size_, max_size_, max_solutions_;
    int sol_count_;
    float accept_score_, max_steps_, min_ss_score_;
    
};

#endif /* defined(__RNAMake__motif_state_search__) */
