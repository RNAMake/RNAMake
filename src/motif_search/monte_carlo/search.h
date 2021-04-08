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
#include <motif_search/solution_topology.h>
#include <motif_search/monte_carlo/scorer.h>

namespace motif_search {
namespace monte_carlo {


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
    apply(
            motif_data_structure::MotifStateGraphOP,
            float) = 0;


    virtual
    float
    score() = 0;

    virtual
    void
    undo(
            motif_data_structure::MotifStateGraphOP) = 0;



public:
    void
    set_temperature(
            float temp) { mc_.set_temperature(temp); }


    void
    scale_temperature(
            float scale) { mc_.scale_temperature(scale); }

protected:
    String name_;
    util::MonteCarlo mc_;
};

typedef std::shared_ptr<Move> MoveOP;
typedef std::vector<MoveOP>   MoveOPs;

class MotifSwapMove : public Move {
public:
    MotifSwapMove(
            ScorerOP scorer,
            SolutionToplogy const & sol_top):
            Move("MotifSwap"),
            scorer_(scorer),
            sol_top_(sol_top) {
        this->mc_ = util::MonteCarlo(1.0f);
        rng_ = util::RandomNumberGenerator();
    }

public:
    bool
    apply(
            motif_data_structure::MotifStateGraphOP msg,
            float current_score) {
        pos_ = rng_.randrange((int)sol_top_.size() - 1);
        new_ms = sol_top_.get_motif_state(pos_);
        last_ms_ = msg->get_node(pos_+1)->data()->cur_state;
        msg->replace_state(pos_+1, new_ms);
        new_score_ = scorer_->score(*msg->last_node()->data()->cur_state->end_states()[1]);
        accept_ = mc_.accept(current_score, new_score_);
        if(accept_) {
            return true;
        }
        else {
            undo(msg);
            return false;
        }

        return true;
    }

    float
    score() { return new_score_; }

    void
    undo(
            motif_data_structure::MotifStateGraphOP msg) {
        msg->replace_state(pos_+1, last_ms_);
    }

private:
    ScorerOP scorer_;
    SolutionToplogy sol_top_;
    float new_score_;
    int accept_, pos_;
    util::RandomNumberGenerator rng_;
    motif::MotifStateOP new_ms, last_ms_;
};


class MoveSet {

};


class Search : public motif_search::Search {
public:
    struct Parameters {
        float accept_score = 10.0f;
        int max_size = 100000;
    };

public:
    Search(
            ScorerOP scorer,
            SolutionToplogy const & sol_top,
            SolutionFilterOP filter):
            motif_search::Search("monte_carlo"),
            scorer_(scorer->clone()),
            sol_top_(sol_top),
            filter_(filter->clone()) {
        finished_ = false;
        setup_options();
        update_var_options();
    }

    ~Search() {}


    motif_search::Search *
    clone() const { return new Search(*this); };

public:

    virtual
    void
    setup(
            ProblemOP p)  {
        msg_ = sol_top_.initialize_solution(p->start);
        scorer_->set_target(p->end, p->target_an_aligned_end);
        lookup_ = p->lookup;
        if(lookup_ != nullptr) {
            using_lookup_ = true;
        }
        stages_ = 50;
        steps_ = 500000;
        stage_ = 0;
    }


    virtual
    void
    //TODO implement or change interface
    start() { }

    virtual
    bool
    //TODO implement or change interface
    finished() {
        return false;
    }

    virtual
    SolutionOP
    next() {
        auto cur_score = scorer_->score(*msg_->last_node()->data()->cur_state->end_states()[1]);
        auto new_score = 0.0;
        auto best_score = cur_score;
        auto mover = std::make_shared<MotifSwapMove>(scorer_, sol_top_);
        auto hot_mover = std::make_shared<MotifSwapMove>(scorer_, sol_top_);
        auto min_mover = std::make_shared<MotifSwapMove>(scorer_, sol_top_);
        hot_mover->set_temperature(100.0f);
        min_mover->set_temperature(0.1f);
        auto temp = 4.5;
        mover->set_temperature(temp);

        // TODO add temperature adjustment to get to 0.235 acceptance
        // TODO add mnimization step
        // TODO add sterics
        auto accept = false;
        auto accepted_steps = 0.0;
        while(stage_ < stages_) {
            accepted_steps = 0.0;
            auto round_best_score = 100000.0f;
            auto best_sol = motif_data_structure::MotifStateGraphOP(nullptr);
            while(step_ < steps_) {
                step_ += 1;
                accept = mover->apply(msg_, cur_score);
                if(!accept) { continue; }
                accepted_steps += 1;
                cur_score = mover->score();
                if(cur_score < parameters_.accept_score) {
                    _get_solution_motif_names(msg_);
                    if(!filter_->accept(motif_names_)) { continue; }
                    auto sol_msg = _get_solution_msg();
                    LOG_DEBUG << "found a solution: " << cur_score;
                    return std::make_shared<Solution>(sol_msg, cur_score);
                }
                if(cur_score < best_score) {
                    best_score = cur_score;
                }
                if(cur_score < round_best_score) {
                    round_best_score = cur_score;
                }

            }
            auto accept_ratio = (float)(accepted_steps / steps_);
            LOG_DEBUG << "stage: " << stage_ << " best_score: " << best_score << " acceptance: " << accept_ratio << " temp: " << temp;
            LOG_DEBUG << "round best score: " << round_best_score;
            // heatup
            for(int i = 0; i < 1000; i++) {
                accept = hot_mover->apply(msg_, cur_score);
                if(accept) {
                    cur_score = hot_mover->score();
                }
            }
            // reset mover
            mover = std::make_shared<MotifSwapMove>(scorer_, sol_top_);
            auto diff = (accept_ratio - 0.235)*10;
            //temp = temp - diff;
            //mover->set_temperature(temp);

            step_ = 0;
            stage_ += 1;
        }

        LOG_DEBUG << "exiting monte carlo search";

        return SolutionOP(nullptr);


    }

private:
    bool
    _steric_clash(
            motif_data_structure::MotifStateGraphOP msg) {
        auto clash = false;
        if (using_lookup_) {
            for (auto const & n : *msg) {
                for(auto const & b : n->data()->cur_state->beads()) {
                    clash = lookup_->clash(b);
                    if (clash) { return true; }
                }
            }
        }
        return false;
    }

    motif_data_structure::MotifStateGraphOP
    _get_solution_msg() {
        auto new_msg = std::make_shared<motif_data_structure::MotifStateGraph>();
        new_msg->set_option_value("sterics", false);
        for(auto const & n : *msg_) {
            if(n->index() == 0) { continue; }
            if(n->index() == 1) {
                new_msg->add_state(std::make_shared<motif::MotifState>(*n->data()->cur_state));
            }
            else {
                new_msg->add_state(std::make_shared<motif::MotifState>(*n->data()->cur_state), -1, n->parent_end_index());
            }

        }
        return new_msg;

    }

    void
    _get_solution_motif_names(
            motif_data_structure::MotifStateGraphOP msg) {
        motif_names_.resize(0);
        for (auto const & n : *msg) {
            motif_names_.push_back(n->data()->name());
        }
    }


protected:

    void
    setup_options() {
        options_.add_option("sterics", true, base::OptionType::BOOL);
        options_.add_option("min_size", 0, base::OptionType::INT);
        options_.add_option("max_size", 1000000, base::OptionType::INT);
        options_.add_option("max_solutions", 1, base::OptionType::INT);
        options_.add_option("accept_score", 10.0f, base::OptionType::FLOAT);
        options_.add_option("return_best", false, base::OptionType::BOOL);
        options_.lock_option_adding();

    }

    void
    update_var_options() {
        parameters_.accept_score = options_.get_float("accept_score");
        parameters_.max_size = options_.get_int("max_size");

    }



private:
    Parameters parameters_;
    ScorerOP scorer_;
    SolutionToplogy sol_top_;
    SolutionFilterOP filter_;
    util::StericLookupNewOP lookup_;
    motif_data_structure::MotifStateGraphOP msg_;
    float score_;
    int stages_, stage_;
    int steps_, step_;
    bool finished_, using_lookup_;
    Strings motif_names_;

};

}
}

#endif //RNAMAKE_NEW_MONTE_CARLO_SEARCH_H
