//
// Created by Joseph Yesselman on 2019-03-30.
//

#include <motif_search/exhaustive/search.h>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// option functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace motif_search::exhaustive {

void
Search::setup_options() {
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("max_node_level", 100, base::OptionType::INT);
    options_.add_option("min_node_level", -1, base::OptionType::INT);
    options_.add_option("min_size", 0, base::OptionType::INT);
    options_.add_option("max_size", 1000000, base::OptionType::INT);
    options_.add_option("max_solutions", 1, base::OptionType::INT);
    options_.add_option("accept_score", 10, base::OptionType::FLOAT);
    options_.add_option("min_ss_score", 10000, base::OptionType::FLOAT);
    options_.add_option("return_best", false, base::OptionType::BOOL);
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
    parameters_.return_best = options_.get_bool("return_best");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// private functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

motif_data_structure::MotifStateGraphOP
Search::_get_graph_from_solution() {
    auto motif_states = enumerator_.all_states();
    auto msg = std::make_shared<motif_data_structure::MotifStateGraph>();
    msg->set_option_value("sterics", false);
    int i = 0, j = 0;
    for(auto const & ms : motif_states) {
        ms->new_uuids();
        if(i == 0) {
            msg->add_state(ms);
        }
        else {
            j = msg->add_state(ms);
            if(j == -1) {
                LOG_ERROR << "something went horribly wrong, cannot build solution";
            }
        }
        i++;
    }
    return msg;
}


}