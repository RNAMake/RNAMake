//
// Created by Joseph Yesselman on 3/17/19.
//

#include "motif_search/path_finding/search.h"

namespace motif_search {
namespace path_finding {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
Search::setup_options() {
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("max_node_level", 100, base::OptionType::INT);
    options_.add_option("min_node_level", 0, base::OptionType::INT);
    options_.add_option("min_size", 0, base::OptionType::INT);
    options_.add_option("max_size", 1000000, base::OptionType::INT);
    options_.add_option("max_solutions", 1, base::OptionType::INT);
    options_.add_option("accept_score", 10, base::OptionType::FLOAT);
    options_.add_option("min_ss_score", 10000, base::OptionType::FLOAT);
    options_.add_option("max_steps", 1000000000, base::OptionType::FLOAT);
    options_.add_option("verbose", true, base::OptionType::BOOL);
    options_.add_option("return_best", false, base::OptionType::BOOL);
    options_.add_option("helix_end", false, base::OptionType::BOOL);
    options_.lock_option_adding();
}

void
Search::update_var_options() {
    parameters_.sterics = options_.get_bool("sterics");
    parameters_.min_size = options_.get_int("min_size");
    parameters_.max_size = options_.get_int("max_size");
    parameters_.max_solutions = options_.get_int("max_solutions");
    parameters_.max_node_level = options_.get_int("max_node_level");
    parameters_.min_node_level = options_.get_int("min_node_level");
    parameters_.accept_score = options_.get_float("accept_score");
    parameters_.min_ss_score = options_.get_float("min_ss_score");
    parameters_.max_steps = options_.get_float("max_steps");
    parameters_.verbose = options_.get_bool("verbose");
    parameters_.helix_end = options_.get_bool("helix_end");
    parameters_.return_best = options_.get_bool("return_best");
}


void
Search::setup(
        motif_search::ProblemOP p)  {

    auto ms = std::make_shared<motif::MotifState>(
            "start", Strings{"start", "start"}, Strings{"", ""},
            structure::BasepairStateOPs {p->start, p->start},
            math::Points(), 0, 0, 0);


}



}
}