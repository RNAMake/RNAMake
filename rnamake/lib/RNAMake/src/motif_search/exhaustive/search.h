//
// Created by Joseph Yesselman on 2019-03-30.
//

#ifndef RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
#define RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H

#include <motif_search/search.h>
#include <motif_search/solution_topology.h>
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
            SolutionToplogy const & sol_top):
            scorer_(scorer->clone()),
            enumerator_(MotifStateEnumerator(sol_top)),
            parameters_(Parameters()),
            options_(base::Options()){
        setup_options();
        update_var_options();
        finished_ = false;
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
        auto count = 0;
        while(!enumerator_.finished()) {
            count += 1;
            //if(count % 1000000 == 0) { std::cout << count << std::endl; }

            current_ = enumerator_.top_state();
            score_ = scorer_->score(*current_->end_states()[1]);

            if(score_ < best) {
                best = score_;
                std::cout << best << " ";
                auto & motif_states = enumerator_.all_states();
                for(auto const & ms : motif_states) {
                    std::cout << ms->name() << " ";
                }
                std::cout << std::endl;
            }

            if(score_ < parameters_.accept_score) {
                auto msg = _get_graph_from_solution();
                return std::make_shared<motif_search::Solution>(msg, score_);
            }
            enumerator_.next();
        }

        finished_ = true;
        return SolutionOP(nullptr);


    }

private:
    motif_data_structure::MotifStateGraphOP
    _get_graph_from_solution();

private:

    void
    setup_options();

    void
    update_var_options();


private:
    Parameters parameters_;
    ScorerOP scorer_;
    MotifStateEnumerator enumerator_;
    base::Options options_;
    motif::MotifStateOP current_;
    float score_;
    bool finished_;

};

}
}



#endif //RNAMAKE_NEW_EXHAUSTIVE_SEARCH_H
