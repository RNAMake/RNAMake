//
//  path_follower.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_search/path_follower.h"

namespace motif_search {

void
PathFollower::set_cmd_options(
        base::CommandLineOptions const & opts) {

    for (auto const & opt : opts) {
        if (!search_.has_option(opt->name())) { continue; }
        if (!opts.is_filled(opt->name())) { continue; }
        if (opt->type() == base::OptionType::INT) {
            search_.set_option_value(opt->name(), opt->get_int());
        } else if (opt->type() == base::OptionType::FLOAT) {
            search_.set_option_value(opt->name(), opt->get_float());
        } else if (opt->type() == base::OptionType::BOOL) {
            search_.set_option_value(opt->name(), opt->get_bool());
        } else if (opt->type() == base::OptionType::STRING) {
            search_.set_option_value(opt->name(), opt->get_string());
        }
    }

    for (auto const & opt: opts) {
        if (!has_option(opt->name())) { continue; }

        if (opt->type() == base::OptionType::INT) {
            set_option_value(opt->name(), opt->get_int());
        } else if (opt->type() == base::OptionType::FLOAT) {
            set_option_value(opt->name(), opt->get_float());
        } else if (opt->type() == base::OptionType::BOOL) {
            set_option_value(opt->name(), opt->get_bool());
        } else if (opt->type() == base::OptionType::STRING) {
            set_option_value(opt->name(), opt->get_string());
        }
    }

}


void
PathFollower::setup_options() {
    options_.add_option("only_one", 0, base::OptionType::INT);
    options_.add_option("max_pathes", 10, base::OptionType::INT);
    options_.add_option("score_diff", .90f, base::OptionType::FLOAT);
    options_.add_option("final_path", false, base::OptionType::BOOL);
    options_.lock_option_adding();
    update_var_options();
}

void
PathFollower::update_var_options() {}


motif_data_structure::MotifTreeOP
PathFollower::next() {
    auto sol = search_.next();
    auto mt = sol->to_motif_tree();

    return mt;

}


MotifStateSearchSolutionOPs
PathFollower::solutions() {
    auto sol = search_.next();
    auto solutions = MotifStateSearchSolutionOPs();

    if (sol == nullptr) { return solutions; }

    search_.setup(start_state_, start_state_);
    search_.set_option_value("max_solutions", options_.get_int("max_pathes"));
    search_.set_option_value("accept_score", sol->score() * 2);
    solutions.push_back(sol);

    while (!search_.finished()) {
        auto sol = search_.next();
        if (sol == nullptr) { continue; }
        if (sol->score() < search_.get_float_option("accept_score")) {
            solutions.push_back(sol);
            search_.set_option_value("accept_score", sol->score());
        }
    }
    return solutions;
}

}
