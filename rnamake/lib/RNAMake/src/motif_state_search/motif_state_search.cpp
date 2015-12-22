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
    options_.add_option(Option("min_node_level", 0));
    options_.add_option(Option("min_size", 0));
    options_.add_option(Option("max_size", 10000));
    options_.add_option(Option("max_solutions", 10));
    options_.add_option(Option("max_steps", LONG_MAX*1.0f));
    options_.add_option(Option("accept_score", 10.0f));
    options_.add_option(Option("min_ss_score", 10000.0f));
    update_var_options();
}

void
MotifStateSearch::update_var_options() {
    sterics_        = options_.option<int>("sterics");
    min_size_       = options_.option<int>("min_size");
    max_size_       = options_.option<int>("max_size");
    max_solutions_  = options_.option<int>("max_solutions");
    max_node_level_ = options_.option<int>("max_node_level");
    min_node_level_ = options_.option<int>("min_node_level");
    accept_score_   = options_.option<float>("accept_score");
    max_steps_      = options_.option<float>("max_steps");
    min_ss_score_   = options_.option<float>("min_ss_score");

}

void
MotifStateSearch::setup(
    const BasepairStateOP & start,
    const BasepairStateOP & end) {
    
    auto start_n = _start_node(start);
    start_n->score(100000);
    test_node_ = std::make_shared<MotifStateSearchNode>(*start_n);
    queue_.push(start_n);
    scorer_->set_target(end);
    no_more_solutions_ = 0;
    sol_count_ = 0;
}

MotifStateSearchSolutionOP
MotifStateSearch::next() {
    auto sol = _search();
    if(sol == nullptr) {
        no_more_solutions_ = 1;
        return nullptr;
    }
    sol_count_++;
    //solutions_.push_back(sol);
    return sol;
    
}

int
MotifStateSearch::finished() {
    if(sol_count_ >= max_solutions_ || no_more_solutions_) {
        return 1;
    }
    else{
        return 0;
    }
}


MotifStateSearchSolutionOP
MotifStateSearch::_search() {
    float score;
    int i = 0;
    int steps = 0;
    int j = 0;
    int pos;
    int path_count = 0;
    int path_size = 0;
    MotifStateSearchNodeOP current, child, current_2;
    int clash;
    float dist;
    float path_score = 0;
    float best_dist;
    Points centers(100);
    int center_count = 0;
    
    while(! queue_.empty() ) {
        current = queue_.top();
        queue_.pop();
        
        
        steps += 1;
        score = scorer_->accept_score(current);
        if(score < accept_score_ && current->ss_score() < min_ss_score_ &&
           current->level() > min_node_level_) {
            auto s = std::make_shared<MotifStateSearchSolution>(current, score);
            return s;
        }
        
        if(current->level()+1 > max_node_level_) { continue; }
        
        possible_children_ = selector_->get_children_ms(current);
        pos = possible_children_.pos();
        test_node_->parent(current);
        Point center;
        i = -1;
        
        for(auto const & end : current->cur_state()->end_states()) {
            i++;
            if(i == 0) { continue; }
            j = -1;
            for(auto const & ms_and_type : possible_children_.motif_states_and_types()) {
                j++;
                
                if(j >= pos) { break; }
                test_node_->replace_ms(ms_and_type.motif_state,
                                       ms_and_type.type);
                
                
                if(current->size() + ms_and_type.motif_state->size() > max_size_) {
                    continue;
                }
                
                aligner_.get_aligned_motif_state(end,
                                                 test_node_->cur_state(),
                                                 test_node_->ref_state());
             
                test_node_->calc_center();

                if(path_.size() > 0) {
                    center_count = 0;
                    current_2 = test_node_;
                    while(current_2 != nullptr) {
                        centers[center_count] = current_2->center();
                        center_count++;
                        current_2 = current_2->parent();
                    }
                    path_count = 0;
                    for(auto const & b : path_) {
                        best_dist = 10000;
                        for(int i = 0; i < center_count; i++) {
                            dist = b.distance(centers[i]);
                            if(dist < best_dist) {
                                best_dist = dist;
                            }
                        }
                        path_score += best_dist;
                    }
                }
                
                //score = scorer_->score(test_node_) + path_score/5;
                //std::cout << path_score << std::endl;
                score = path_score;
                if(score > current->score()) { continue; }
                
                if(sterics_) {
                    clash = 0;
                    for(auto const & b1 : beads_) {
                        for(auto const & b2 : test_node_->cur_state()->beads()) {
                            dist = b1.distance(b2);
                            if(dist < 2.5) { clash = 1; break; }
                        }
                    }
                    if(clash) { continue; }
                    
                    current_2 = test_node_->parent();
                    while (current_2 != nullptr) {
                        for(auto const & b1 : test_node_->cur_state()->beads()) {
                            for(auto const & b2 : current_2->cur_state()->beads()) {
                                dist = b1.distance(b2);
                                if(dist < 2.5) { clash = 1; break; }
                            }
                            if(clash) { break; }
                        }
                        current_2 = current_2->parent();
                    }
                }
                
             
                
                child = std::make_shared<MotifStateSearchNode>(*test_node_);
                child->update();
                child->score(score);
                queue_.push(child);
            }
        }
        
    }
    
    return nullptr;

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