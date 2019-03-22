//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SEARCH_H
#define RNAMAKE_NEW_PATH_FINDING_SEARCH_H

#include <base/option.h>
#include <motif/motif_state_aligner.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_search/search.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/selector.h>
#include <motif_search/path_finding/node.h>

namespace motif_search {
namespace path_finding {

struct Solution : public motif_search::Solution {
    inline
    Solution(
            motif_data_structure::MotifStateGraphOP n_graph,
            float n_score):
            graph(n_graph),
            score(n_score) {}

    motif_data_structure::MotifStateGraphOP graph;
    float score;
};

typedef std::shared_ptr<Solution> SolutionOP;

struct Parameters {
    bool sterics, verbose, helix_end;
    int max_node_level, min_size, max_size, max_solutions;
    int sol_count, min_node_level;
    float accept_score, min_ss_score, max_steps;
    bool return_best;
};

class Search : public motif_search::Search<Solution> {
public:

    Search(
            ScorerOP scorer,
            SelectorOP selector):
            motif_search::Search<Solution>(),
            scorer_(scorer->clone()),
            selector_(selector->clone()),
            aligner_(motif::MotifStateAligner()) {
        parameters_ = Parameters();
        setup_options();
        update_var_options();
    }

public:
    void
    setup(
           motif_search::ProblemOP);

    void
    start() {}

    void
    finished() {}

    SolutionOP
    next();


public: //option wrappers

    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }

    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }

    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }

    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }

    template<typename T>
    void
    set_option_value(
            String const & name,
            T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }

protected:

    void
    update_var_options();

    void
    setup_options();

private:

    motif_data_structure::MotifStateGraphOP
    _graph_from_node(
            NodeOP);

    inline
    bool
    _accept_node(
            Node const & n) {
        if(n.ss_score() > parameters_.min_ss_score) { return false; }
        if(n.level() < parameters_.min_node_level)  { return false; }
        // this is bad ... at motif_type to MotifState? -- JDY
        if(parameters_.helix_end && n.state()->name()[0] != 'H') { return false;}
        return true;
    }

    bool
    _steric_clash(
            motif::MotifState const &,
            Node const &);


private:
    ScorerOP scorer_;
    SelectorOP selector_;
    Parameters parameters_;
    NodeQueue queue_;
    motif::MotifStateAligner aligner_;
    util::StericLookupNewOP lookup_;

    base::Options options_;
    bool using_lookup_;
    bool enumerating_;
};

}
}

#endif //RNAMAKE_NEW_PATH_FINDING_SEARCH_H
