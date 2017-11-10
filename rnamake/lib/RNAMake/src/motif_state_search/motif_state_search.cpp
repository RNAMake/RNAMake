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
    options_.add_option("sterics", true, OptionType::BOOL);
    options_.add_option("max_node_level", 100, OptionType::INT);
    options_.add_option("min_node_level", 0, OptionType::INT);
    options_.add_option("min_size", 0, OptionType::INT);
    options_.add_option("max_size", 1000000, OptionType::INT);
    options_.add_option("max_solutions", 1, OptionType::INT);
    options_.add_option("accept_score", 10, OptionType::FLOAT);
    options_.add_option("min_ss_score", 10000, OptionType::FLOAT);
    options_.add_option("max_steps", 1000000000, OptionType::FLOAT);
    options_.add_option("verbose", true, OptionType::BOOL);
    options_.add_option("return_best", false, OptionType::BOOL);
    options_.add_option("helix_end", false, OptionType::BOOL);
    
    //for making a movie
    options_.add_option("save_midpoints", false, OptionType::BOOL);
    options_.add_option("save_midpoints_file", "midpoints.dat", OptionType::STRING);
    options_.add_option("save_midpoints_freq", 100, OptionType::INT);
    
    options_.lock_option_adding();
    
    update_var_options();
}

void
MotifStateSearch::update_var_options() {
    sterics_        = options_.get_bool("sterics");
    min_size_       = options_.get_int("min_size");
    max_size_       = options_.get_int("max_size");
    max_solutions_  = options_.get_int("max_solutions");
    max_node_level_ = options_.get_int("max_node_level");
    min_node_level_ = options_.get_int("min_node_level");
    accept_score_   = options_.get_float("accept_score");
    min_ss_score_   = options_.get_float("min_ss_score");
    max_steps_      = options_.get_float("max_steps");
    verbose_        = options_.get_bool("verbose");
    helix_end_      = options_.get_bool("helix_end");
    
}

void
MotifStateSearch::setup(
    BasepairStateOP const & start,
    BasepairStateOP const & end) {
    
    auto start_n = _start_node(start);
    start_n->score(1000000000);
    queue_ =  MotifStateSearchNodeQueue();
    test_node_ = std::make_shared<MotifStateSearchNode>(*start_n);
    queue_.push(start_n);
    scorer_->set_target(end);
    no_more_solutions_ = 0;
    sol_count_ = 0;
}

void
MotifStateSearch::reset() {
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
    return sol;
    
}

int
MotifStateSearch::finished() {
    if(sol_count_ >= max_solutions_ || no_more_solutions_) {
        if(verbose_) {
            std::cout << "maxed out solutions " << sol_count_ << " " << max_solutions_ <<  std::endl;
        }
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
    MotifStateSearchNodeOP current, child, current_2;
    int clash;
    float dist = 1000;
    float best = 1000000000;
    MotifStateSearchSolutionOP best_sol;
    
    std::ofstream midpoint_file;
    auto midpoints = get_bool_option("save_midpoints");
    int midpoint_freq = get_int_option("save_midpoints_freq");
    
    if(midpoints) { midpoint_file.open(get_string_option("save_midpoints_file")); }
    
    while(! queue_.empty() ) {
        current = queue_.top();
        queue_.pop();
        
        
        steps += 1;
        
        if(steps > max_steps_) {
            if(verbose_) {
                std::cout << "MOTIF STATE SEARCH: reached max steps" << std::endl;
            }
            return best_sol;
        }
        
        if(midpoints && steps % midpoint_freq == 0) {
            auto sol = MotifStateSearchSolution(current, score);
            midpoint_file << sol.to_motif_tree()->topology_to_str() << std::endl;
            
        }

        score = scorer_->accept_score(current);

        if(score < best) {
            best = score;
            best_sol = std::make_shared<MotifStateSearchSolution>(current, score);

            if(verbose_) {
                std::cout << "MOTIF STATE SEARCH: best_score=" << best << " motifs_in_solution=";
                std::cout << current->level()-1 << " steps=" << steps << std::endl;
                
                //std::cout << best << " " << accept_score_ << " " << current->level() << " " << steps << " " << max_steps_ << std::endl;
            }
        }
        
        if(score < accept_score_ && current->ss_score() < min_ss_score_ &&
           current->level() > min_node_level_) {
            if(helix_end_ && current->ref_state()->name()[0] != 'H') {
                continue;
            }


            if(verbose_) {
                std::cout << "MOTIF STATE SEARCH: found a solution!" << std::endl;
            }
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
                //std::cout << current->level() << " " << current->ntype() << " " << ms_and_type.type << std::endl;
                test_node_->replace_ms(ms_and_type.motif_state,
                                       ms_and_type.type);
                
                
                if(current->size() + ms_and_type.motif_state->size() > max_size_) {
                    continue;
                }
                
                aligner_.get_aligned_motif_state(end,
                                                 test_node_->cur_state(),
                                                 test_node_->ref_state());
                test_node_->update(i);
                
                score = scorer_->score(test_node_);
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
                    
                    if(using_lookup_) {
                        for(auto const & p : test_node_->cur_state()->beads()) {
                            //std::cout << p << " " << ms_and_type.motif_state->name() <<  std::endl;
                            clash = lookup_.clash(p);
                            if(clash) { break; }
                        }
                        if(clash) { continue; }
                    }
                    
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
                    if(clash) { continue; }

                }
                
             
                
                child = std::make_shared<MotifStateSearchNode>(*test_node_);
                child->update(i);
                child->score(score);
                queue_.push(child);
            }
        }
        
    }
    if(verbose_) {
        std::cout << "MOTIF STATE SEARCH: ran out of options" << std::endl;
    }
    
    if(midpoints) { midpoint_file.close(); }

    
    return best_sol;

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
