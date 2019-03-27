//
// Created by Joseph Yesselman on 3/23/19.
//

#ifndef RNAMAKE_NEW_MONTE_CARLO_SEARCH_H
#define RNAMAKE_NEW_MONTE_CARLO_SEARCH_H

#include <base/option.h>
#include <util/monte_carlo.h>
#include <motif/motif_state_aligner.h>
#include <motif_data_structure/motif_state_graph.hpp>
#include <motif_search/search.h>
#include <motif_search/solution_filter.h>
#include <motif_search/monte_carlo/scorer.h>
#include <motif_search/monte_carlo/solution_topology.h>

namespace motif_search {
namespace monte_carlo {

struct State {
public:
    inline
    State(
            motif_data_structure::MotifStateGraph const & n_msg,
            SolutionToplogy const & n_sol_topology,
            util::StericLookupNewOP n_lookup,
            ScorerOP n_scorer):
            msg(n_msg),
            sol_topology(n_sol_topology),
            lookup(n_lookup),
            scorer(n_scorer->clone()) {
        if(lookup->size() == 0) { using_lookup = false; }
        else                     { using_lookup = true; }
        score = get_score();
    }

public:
    inline
    float
    get_score() {
        auto best_score = 1000.0;
        auto cur_score = 0.0;
        for(auto const & nie : sol_topology.get_solution_nie()) {
            cur_score = scorer->score(*msg.get_node(nie.node_index)->data()->get_end_state(nie.edge_index));
            if(cur_score < best_score) {
                best_score = cur_score;
            }
        }
        return best_score;
    }


    bool
    steric_clash() {
        auto clash = false;
        if(using_lookup) {
            for(auto const & n : msg) {
                clash = lookup->clash(n->data()->cur_state->beads());
                if(clash) { return true; }
            }
        }
        return false;
    }

public:
    motif_data_structure::MotifStateGraph msg;
    SolutionToplogy sol_topology;
    util::StericLookupNewOP lookup;
    ScorerOP scorer;
    float score;
    bool using_lookup;
};

typedef std::shared_ptr<State> StateOP;

class Move {
public:
    Move(
            String const & name):
            name_(name),
            mc_(util::MonteCarlo(3.0)) {}

    virtual
    ~Move() {}

public:
    virtual
    bool
    attempt(
            State &) = 0;

protected:
    bool clash_;
    String name_;
    util::MonteCarlo mc_;
};

class MotifSwapMove : Move {

    MotifSwapMove() : Move("MotifSwap") {
        this->mc_ = util::MonteCarlo(10.0f);
    }

private:
    float new_score_;
    int accept_, pos_;
    motif::MotifStateOP last_ms_;
};




class Search : public motif_search::Search {
public:
    struct Parameters {

    };

public:
    Search(
            ScorerOP scorer,
            SolutionToplogy const & sol_top):
            scorer_(scorer->clone()),
            sol_top_(sol_top) {

    }

    ~Search() {}


    motif_search::Search *
    clone() const { return new Search(*this); };

public:

    virtual
    void
    setup(
            ProblemOP p)  {
        auto msg = sol_top_.initialize_solution(p->start);
        //state_ = std::make_shared<State>();

    }


    virtual
    void
    start() {}

    virtual
    void
    finished() {}

    virtual
    SolutionOP
    next() { return SolutionOP(nullptr); }


private:

private:
    StateOP state_;
    ScorerOP scorer_;
    SolutionToplogy sol_top_;
    base::Options options_;

};

}
}

#endif //RNAMAKE_NEW_MONTE_CARLO_SEARCH_H
