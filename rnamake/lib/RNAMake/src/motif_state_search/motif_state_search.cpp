//
//  motif_state_search.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <climits>

#include "motif/motif_state.h"
#include "motif_state_search/motif_state_search.h"

void
MotifStateSearch::setup_options() {
    options_ = Options();
    options_.add_option(Option("sterics", 1));
    options_.add_option(Option("max_node_level", 15));
    options_.add_option(Option("min_size", 0));
    options_.add_option(Option("max_size", 10000));
    options_.add_option(Option("max_solutions", 10));
    options_.add_option(Option("max_steps", LONG_MAX*1.0f));
    options_.add_option(Option("accept_score", 10.0f));
    update_var_options();
}

void
MotifStateSearch::update_var_options() {
    sterics_        = options_.option<int>("sterics");
    min_size_       = options_.option<int>("min_size");
    max_size_       = options_.option<int>("max_size");
    max_solutions_  = options_.option<int>("max_solutions");
    max_node_level_ = options_.option<int>("max_node_level");
    accept_score_   = options_.option<float>("accept_score");
    max_steps_      = options_.option<float>("max_steps");

}

MotifStateSearchSolutionOPs
MotifStateSearch::search(
    BasepairStateOP const & start,
    BasepairStateOP const & end) {
    
    auto start_n = _start_node(start);
    auto test_node = std::make_shared<MotifStateSearchNode>(start_n->copy());
    MotifStateSearchNodeOP current, child;
    queue_.push(start_n);
    scorer_->set_target(end);
    
    float score;
    int i = 0;
    int steps = 0;
    int j = 0;
    int pos;
    
    while(! queue_.empty() ) {
        current = queue_.top();
        queue_.pop();
        
        if(steps % 1000 == 0) {
        //std::cout << steps << " " << current->level() << " " << current->score() << std::endl;
        }
        
        steps += 1;
        score = scorer_->accept_score(current);
        if(score < accept_score_) {
            auto s = std::make_shared<MotifStateSearchSolution>(current, score);
            solutions_.push_back(s);
            if(solutions_.size() >= max_solutions_) { return solutions_; }
        }
        
        if(current->level()+1 > max_node_level_) { continue; }
        
        possible_children_ = selector_->get_children_ms(current);
        pos = possible_children_.pos();
        test_node->parent(current);
        i = -1;
                
        for(auto const & end : current->cur_state()->end_states()) {
            i++;
            if(i == 0) { continue; }
            j = -1;
            for(auto const & ms_and_type : possible_children_.motif_states_and_types()) {
                j++;
                
                if(j >= pos) { break; }
                test_node->replace_ms(ms_and_type.motif_state,
                                      ms_and_type.type);
                
                aligner_.get_aligned_motif_state(end,
                                                 test_node->cur_state(),
                                                 test_node->ref_state());
                
                score = scorer_->score(test_node);
                if(score > current->score()) { continue; }
                child = std::make_shared<MotifStateSearchNode>(test_node->copy());
                child->update();
                child->score(score);
                queue_.push(child);
            }
        }
        
    }
    
    return solutions_;
    
}



MotifStateSearchNodeOP
MotifStateSearch::_start_node(
    BasepairStateOP const & start_bp) {
    
    auto ms = std::make_shared<MotifState>("start", Strings {"start", "start"},
                                           Strings {"", ""},
                                           BasepairStateOPs { start_bp, start_bp},
                                           Points(), 0, 0, 0);
    
    auto n = std::make_shared<MotifStateSearchNode>(ms, nullptr, -1, -1);
    n->setup_node_type_usage(selector_->size());
    
    return n;
}