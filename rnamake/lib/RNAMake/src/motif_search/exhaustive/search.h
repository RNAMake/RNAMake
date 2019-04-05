//
// Created by Joseph Yesselman on 2019-03-30.
//

#ifndef RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
#define RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H

#include <motif_search/search.h>
#include <motif_search/solution_topology.h>
#include <motif_search/solution_filter.h>
#include <motif_search/exhaustive/scorer.h>
#include <motif_search/exhaustive/motif_state_enumerator.h>

namespace motif_search {
namespace exhaustive {

class Search : public motif_search::Search {
public:
    struct Parameters {
        bool sterics, helix_end;
        int max_node_level, min_size, max_size, max_solutions, min_node_level;
        float accept_score, min_ss_score, max_steps;
        bool return_best;
    };

public:
    Search(
            ScorerOP scorer,
            SolutionToplogy const & sol_top,
            SolutionFilterOP filter):
            motif_search::Search("exhaustive"),
            scorer_(scorer->clone()),
            enumerator_(MotifStateEnumerator(sol_top)),
            filter_(filter->clone()),
            parameters_(Parameters()){
        finished_ = false;
        setup_options();
        update_var_options();
    }

    ~Search() {}

    motif_search::Search *
    clone() const { return new Search(*this); };

public:

    void
    setup(
            ProblemOP p)  {
        enumerator_.start(p->start);
        scorer_->set_target(p->end, p->target_an_aligned_end);
    }

    void
    start() {}

    bool
    finished() {
        return finished_;
    }

    SolutionOP
    next() {
        auto best = 1000.0f;
        while(!enumerator_.finished()) {
            current_ = enumerator_.top_state();
            score_ = scorer_->score(*current_->end_states()[1]);

            if(score_ < parameters_.accept_score) {
                auto msg = _get_graph_from_solution();
                enumerator_.next();
                _get_solution_motif_names(msg);
                if(!filter_->accept(motif_names_)) { continue; }

                return std::make_shared<motif_search::Solution>(msg, score_);
            }
            enumerator_.next();
        }

        finished_ = true;
        return SolutionOP(nullptr);


    }

private:
    void
    _get_solution_motif_names(
            motif_data_structure::MotifStateGraphOP msg) {
        motif_names_.resize(0);
        for (auto const & n : *msg) {
            motif_names_.push_back(n->data()->name());
        }
    }

    motif_data_structure::MotifStateGraphOP
    _get_graph_from_solution();

protected:

    void
    setup_options();

    void
    update_var_options();


private:
    Parameters parameters_;
    ScorerOP scorer_;
    SolutionFilterOP filter_;
    MotifStateEnumerator enumerator_;
    motif::MotifStateOP current_;
    float score_;
    bool finished_;
    Strings motif_names_;

};

}
}



#endif //RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
