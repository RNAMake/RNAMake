//
//  path_follower.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_state_search/path_follower.h"

void
PathFollower::set_cmd_options(
    CommandLineOptions const & opts) {
    
    for(auto const & opt : opts) {
        if(! search_.has_option(opt->name())) { continue; }
        if(! opts.is_filled(opt->name())) { continue; }
        if     (opt->type() == OptionType::INT) {
            search_.set_option_value(opt->name(), opt->get_int());
        }
        else if(opt->type() == OptionType::FLOAT) {
            search_.set_option_value(opt->name(), opt->get_float());
        }
        else if(opt->type() == OptionType::BOOL) {
            search_.set_option_value(opt->name(), opt->get_bool());
        }
        else if(opt->type() == OptionType::STRING) {
            search_.set_option_value(opt->name(), opt->get_string());
        }
    }

    for(auto const & opt: opts) {
        if(! has_option(opt->name())) { continue; }
        
        if     (opt->type() == OptionType::INT) {
            set_option_value(opt->name(), opt->get_int());
        }
        else if(opt->type() == OptionType::FLOAT) {
            set_option_value(opt->name(), opt->get_float());
        }
        else if(opt->type() == OptionType::BOOL) {
            set_option_value(opt->name(), opt->get_bool());
        }
        else if(opt->type() == OptionType::STRING) {
            set_option_value(opt->name(), opt->get_string());
        }
    }
    
}


void
PathFollower::setup_options() {
    options_.add_option("only_one", 0, OptionType::INT);
    options_.add_option("max_pathes", 10, OptionType::INT);
    options_.add_option("score_diff", .95f, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
PathFollower::update_var_options() {}


MotifTreeOP
PathFollower::next() {
    auto sol = search_.next();
    auto mt = sol->to_motif_tree();
    
    return mt;
    
}


MotifStateSearchSolutionOPs
PathFollower::solutions() {
    auto sol = search_.next();
    auto solutions = MotifStateSearchSolutionOPs();
    
    search_.setup(start_->state(), start_->state());
    search_.set_option_value("max_solutions", options_.get_int("max_pathes"));
    search_.set_option_value("accept_score", sol->score()*options_.get_float("score_diff"));

    while(!search_.finished()) {
        auto sol = search_.next();
        if(sol->score() < search_.get_float_option("accept_score")) {
            solutions.push_back(sol);
        }
    }
    
    return solutions;
    
}
