//
//  motif_tree_state_search.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_search.h"
#include "motif_tree_state_search_scorer.h"
#include "motif_tree_state.h"

/*
 MotifTreeStateSearchOptions functions
 **********************************************************************************
 */

void
MotifTreeStateSearchOptions::set_defaults() {
    numeric_options_["accept_score"   ] = 10.0;
    numeric_options_["max_steps"      ] = 10000000000;
    numeric_options_["max_node_level" ] = 20;
    numeric_options_["sterics"        ] = 1;
    numeric_options_["max_n_solutions"] = 10;
    numeric_options_["accept_ss_score"] = 1000;
    numeric_options_["max_size"       ] = 1000;
    numeric_options_["scorer"         ] = 0;
}

void
MotifTreeStateSearchOptions::numeric(
    String const option,
    float const value) {
    
    if(numeric_options_.find(option) == numeric_options_.end()) {
        std::cout << "cannot find option " << option << " in MotifTreeStateSearchOptions::numeric" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    numeric_options_[option] = value;
}

void
MotifTreeStateSearchOptions::numerics(
    StringFloatMap const & options) {
    
    for(auto const & kv : options) {
        numeric(kv.first,kv.second);
    }
    
}

const
float
MotifTreeStateSearchOptions::numeric(
    String const option) {
    
    if(numeric_options_.find(option) == numeric_options_.end()) {
        std::cout << "cannot find option " << option << " in MotifTreeStateSearchOptions" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return numeric_options_[option];
}


MotifTreeStateSearchSolutionOPs const &
MotifTreeStateSearch::search(
    BasepairStateOP const & start,
    BasepairStateOP const & end,
    MotifTreeStateSelectorOP const & node_selector) {
    
    if (node_selector != NULL) { node_selector_ = node_selector; }
    scorer_ = MotifTreeStateSearchScorerOP(new MTSS_GreedyBestFirstSearch ( end ));
    MotifTreeStateSearchNodeOP start_node = _get_start_node(start);
    
    MotifTreeStateSearchNodeOP test_node (new MotifTreeStateSearchNode(start_node));
    MotifTreeStateSearchNodeOP current = NULL;
    BasepairStateOPs parent_ends;
    MTSTypePairs mts_type_pairs;
    queue_.push(start_node);
    
    //get user defined options
    //float max_steps       = options_.numeric("max_steps");
    float max_size        = options_.numeric("max_size");
    float min_score       = options_.numeric("accept_score");
    float max_node_level  = options_.numeric("max_node_level");
    int   sterics         = options_.numeric("sterics");
    int   max_n_solutions = options_.numeric("max_n_solutions");
    float score = 0.0f;
    int fail = 0;

    while (!queue_.empty() ) {
        current = queue_.top();
        queue_.pop();
        
        steps_++;
        
        if(current->size() > max_size)         { continue; }
        score = scorer_->accept_score(current);
        if ( score < min_score) {
            fail = 0;
            if(!fail) {
                MotifTreeStateSearchSolutionOP solution ( new MotifTreeStateSearchSolution(current, score) );
                solutions_.push_back(solution);
                if(solutions_.size() == max_n_solutions) { return solutions_; }
            }
        }
        
        if(current->level() == max_node_level) { continue; }
        test_node->parent(current);
        test_node->level(current->level() + 1);
        parent_ends = current->active_states();
        mts_type_pairs = node_selector_->get_children_mts(current);
        
        for(auto const & mts_type_pair : mts_type_pairs) {
            for (auto const & end : parent_ends) {
                test_node->replace_mts(mts_type_pair.first);
                aligner_.transform_state(end, current, test_node);
                score = scorer_->score(test_node);
                if (score > current->score()) { continue; }
                if (sterics) {
                    aligner_.transform_beads(test_node);
                    if(test_node->steric_clash()) { continue; }
                }
                test_node->lib_type(mts_type_pair.second);
                test_node->score(score);
                test_node->update_node_count();
                queue_.push(MotifTreeStateSearchNodeOP(new MotifTreeStateSearchNode(test_node)));
            }
        }
        
        //if(steps_ > 10) { break; }
    }
    
    return solutions_;
}


MotifTreeStateSearchNodeOP
MotifTreeStateSearch::_get_start_node(
    BasepairStateOP const & start) {
    
    BasepairStateOPs states(1);
    states[0] = start;
    MotifTreeStateOP mts (new MotifTreeState("start", 1, 0, 0, Points(), states, 0, ""));
    MotifTreeStateSearchNodeOP start_node ( new MotifTreeStateSearchNode(mts, NULL, -1));
    //start_node->setup_node_count((int)node_selector_->nodes().size());
    start_node->setup_node_count(1);

    return start_node;
}







