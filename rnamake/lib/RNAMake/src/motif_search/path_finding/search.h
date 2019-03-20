//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SEARCH_H
#define RNAMAKE_NEW_PATH_FINDING_SEARCH_H

#include <base/option.h>
#include <motif_search/search.h>
#include <motif_search/path_finding/scorer.h>
#include <motif_search/path_finding/selector.h>
#include <motif_search/path_finding/node.h>

namespace motif_search {
namespace path_finding {

struct Solution : public motif_search::Solution {

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
            Scorer const & scorer,
            Selector const & selector):
            motif_search::Search<Solution>(),
            scorer_(scorer),
            selector_(selector) {

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
    next() {
        return SolutionOP(nullptr);
    }


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

private:
    Scorer const & scorer_;
    Selector const & selector_;
    Parameters parameters_;
    NodeQueue queue_;


    base::Options options_;
    int using_lookup_;
    bool enumerating_;
};

}
}

#endif //RNAMAKE_NEW_PATH_FINDING_SEARCH_H
