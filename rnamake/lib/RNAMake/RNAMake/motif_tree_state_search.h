//
//  motif_tree_state_search.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_search__
#define __RNAMake__motif_tree_state_search__

#include <stdio.h>
#include <queue>
#include "motif_type.h"
#include "motif_tree_state_node_aligner.h"
#include "motif_tree_state_search_node.h"
#include "motif_tree_state_search_solution.h"
#include "motif_tree_state_selector.h"
#include "motif_tree_state_search_scorer.h"

class MotifTreeStateSearchNodeOPCompare {
public:
    bool
    operator () (
                 MotifTreeStateSearchNodeOP const & node1,
                 MotifTreeStateSearchNodeOP const & node2) {
        
        if (node1->score() > node2->score()) { return true;  }
        else     							 { return false; }
    }
};

typedef std::priority_queue<MotifTreeStateSearchNodeOP, MotifTreeStateSearchNodeOPs, MotifTreeStateSearchNodeOPCompare> MotifTreeStateSearchNodeQueue;

class MotifTreeStateSearchOptions{
public:
    MotifTreeStateSearchOptions():
    numeric_options_ ( StringFloatMap() ) {
        set_defaults();
    }
    
    ~MotifTreeStateSearchOptions()
    {}
    
public:
    void
    numerics(
             StringFloatMap const &);
    
public:
    void
    numeric(
            std::string const,
            float const);
    
    const
    float
    numeric(
            std::string const);
    
    const
    StringFloatMap
    numeric_options() { return numeric_options_; }
    
private:
    void
    set_defaults();
    
private:
    StringFloatMap numeric_options_;
};

class MotifTreeStateSearch {
public:
    MotifTreeStateSearch():
    queue_ ( MotifTreeStateSearchNodeQueue() ),
    solutions_ ( MotifTreeStateSearchSolutionOPs() ),
    aligner_ ( MotifTreeStateNodeAligner() ),
    options_ ( MotifTreeStateSearchOptions() ),
    steps_ ( 0 ),
    base_beads_ ( Points() )
    {
        node_selector_ = MotifTreeStateSelectorOP ( new MotifTreeStateSelector ( default_selector(MotifTypes())) );
        scorer_ = MotifTreeStateSearchScorerOP(new MTSS_GreedyBestFirstSearch ());

    }
    
    ~MotifTreeStateSearch() {}
    
public:
    
    MotifTreeStateSearchSolutionOPs const &
    search(
        BasepairStateOP const &,
        BasepairStateOP const &,
        MotifTreeStateSelectorOP const & node_selector = NULL,
        MotifTreeStateSearchScorerOP const & scorer = NULL);
    
    void
    reset() {
        queue_ = MotifTreeStateSearchNodeQueue();
        steps_ = 0;
        solutions_.resize(0);
    }
    
public: //setters
    
    inline
    void
    set_numeric_option(
        String const option,
        float const value) {
        options_.numeric(option, value);
    }
    
    inline
    void
    set_numeric_options(
        StringFloatMap const & noptions) {
        options_.numerics(noptions);
    }
    
    inline
    void
    base_beads(Points const & nbase_beads) { base_beads_ = nbase_beads; }
    
private:
    
    MotifTreeStateSearchNodeOP
    _get_start_node(
        BasepairStateOP const &);
    
private:
    MotifTreeStateSearchNodeQueue queue_;
    MotifTreeStateSearchSolutionOPs solutions_;
    MotifTreeStateSearchScorerOP scorer_;
    MotifTreeStateSelectorOP node_selector_;
    MotifTreeStateNodeAligner aligner_;
    MotifTreeStateSearchOptions options_;
    Points base_beads_;
    
    int steps_;
};


#endif /* defined(__RNAMake__motif_tree_state_search__) */
