//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SEARCH_H
#define RNAMAKE_NEW_PATH_FINDING_SEARCH_H

#include <base/option.h>
#include <motif/motif_state_aligner.h>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/selector.h>
#include <motif_search/path_finding/node.h>

namespace motif_search {
namespace path_finding {

struct Parameters {
    bool sterics, helix_end;
    int max_node_level, min_size, max_size, max_solutions, min_node_level;
    float accept_score, min_ss_score, max_steps;
    bool return_best;
};

class Search : public motif_search::Search {
public:

    Search(
            ScorerOP scorer,
            SelectorOP selector,
            SolutionFilterOP solution_filter):
            motif_search::Search(),
            scorer_(scorer->clone()),
            selector_(selector->clone()),
            solution_filter_(solution_filter->clone()),
            aligner_(motif::MotifStateAligner()) {
        motif_names_ = Strings();
        motif_names_.reserve(100);
        parameters_ = Parameters();
        setup_options();
        update_var_options();
    }

    ~Search() {}

    motif_search::Search *
    clone() const { return new Search(*this); };

public:
    void
    setup(
           motif_search::ProblemOP);

    void
    start() {}

    void
    finished() {}

    motif_search::SolutionOP
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

    void
    _get_solution_motif_names(
            NodeOP,
            Strings &);


private:
    ScorerOP scorer_;
    SelectorOP selector_;
    SolutionFilterOP solution_filter_;
    Parameters parameters_;
    NodeQueue queue_;
    motif::MotifStateAligner aligner_;
    util::StericLookupNewOP lookup_;

    base::Options options_;
    Strings motif_names_;
    bool using_lookup_, enumerating_;
};

}
}

#endif //RNAMAKE_NEW_PATH_FINDING_SEARCH_H
