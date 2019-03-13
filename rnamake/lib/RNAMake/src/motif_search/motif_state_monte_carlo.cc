//
// Created by Joseph Yesselman on 11/4/17.
//

#include <fstream>
#include <motif_search/motif_state_monte_carlo.h>

namespace motif_search {

void
MotifStateMonteCarlo::setup_options() {
    options_.add_option("sterics", false, base::OptionType::BOOL);
    options_.add_option("max_solutions", 1, base::OptionType::INT);
    options_.add_option("accept_score", 5, base::OptionType::FLOAT);
    options_.add_option("stages", 100, base::OptionType::INT);
    options_.add_option("steps", 500000, base::OptionType::INT);
    options_.lock_option_adding();

    update_var_options();
}

void
MotifStateMonteCarlo::update_var_options() {
    sterics_ = options_.get_bool("sterics");
    max_solutions_ = options_.get_int("max_solutions");
    accept_score_ = options_.get_float("accept_score");
    stages_ = options_.get_int("stages");
    steps_ = options_.get_int("steps");

}

void
MotifStateMonteCarlo::setup(
        motif_data_structure::MotifStateGraphOP msg,
        int ni,
        int nj,
        int ei,
        int ej,
        bool target_an_aligned_end) {

    msg_ = msg;
    msg_->set_option_value("sterics", false);
    ni_ = ni;
    nj_ = nj;
    ei_ = ei;
    ej_ = ej;
    target_an_aligned_end_ = target_an_aligned_end;

    end_ = msg_->get_node(nj_)->data()->get_end_state(ej_);
    end_flip_ = structure::BasepairStateOP(new structure::BasepairState(end_->copy()));
    end_flip_->flip();

    org_num_ = msg_->size();

    auto i = -1;
    for (auto const & motif_states : mses_) {
        i++;
        auto ms = motif_states[rng_.randrange((int) motif_states.size() - 1)];
        if (i == 0) {
            msg_->add_state(ms, ni_, ei_);
        } else {
            msg_->add_state(ms);
        }
    }

    stage_ = 0;
    step_ = 0;
    enumerating_ = false;
    seen_ = StringIntMap();
}

MotifStateMonteCarloSolutionOP
MotifStateMonteCarlo::next() {
    // do not call run() now, getting one solution at a time
    enumerating_ = true;

    auto score_num = 0;
    auto pos = 0;
    auto cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
    auto new_score = 0.0;
    auto best_score = cur_score;
    mc_.set_temperature(1.0);
    while (stage_ < stages_) {
        while (step_ < steps_) {
            step_ += 1;
            new_score = perform_motif_swap(cur_score);
            // monte carlo move accepted
            if (new_score > 0) { cur_score = new_score; }
            else { continue; }

            if (new_score < best_score) { best_score = new_score; }

            // accept solution?
            if (new_score < accept_score_) {
                if (_seen_solution(msg_)) { continue; }
                auto mg = msg_->to_motif_graph();
                auto last_n = msg_->last_node()->index();
                auto nj_name = msg_->get_node(nj_)->data()->end_name(ej_);
                auto last_name = msg_->get_node(last_n)->data()->end_name(1);
                mg->add_connection(nj_, last_n, nj_name, last_name);
                return std::make_shared<MotifStateMonteCarloSolution>(mg, new_score);
            }

        }
        std::cout << "stage: " << stage_ << " best_score: " << best_score << " cur_score: " << cur_score << std::endl;

        // heat back up
        /*mc_.set_temperature(1000.0);
        for(int i = 0; i < 100; i++) {
            new_score = perform_motif_swap(cur_score);
            // monte carlo move accepted
            if(new_score > 0) { cur_score = new_score; }
            else              { continue; }
        }
        mc_.set_temperature(1.0);
        */
        int j = 0;
        for (int i = org_num_; i < org_num_ + mses_.size(); i++) {
            auto new_ms = mses_[j][rng_.randrange((int) mses_[pos].size() - 1)];
            msg_->replace_state(i, new_ms);
            j++;
        }

        cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);

        step_ = 0;
        stage_ += 1;
    }

    return MotifStateMonteCarloSolutionOP(nullptr);
}


MotifStateMonteCarloSolutionNewOP
MotifStateMonteCarlo::next_state() {
    // do not call run() now, getting one solution at a time
    enumerating_ = true;

    auto score_num = 0;
    auto pos = 0;
    auto cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
    auto new_score = 0.0;
    auto best_score = cur_score;
    mc_.set_temperature(1.0);
    while (stage_ < stages_) {
        while (step_ < steps_) {
            step_ += 1;
            new_score = perform_motif_swap(cur_score);
            // monte carlo move accepted
            if (new_score > 0) { cur_score = new_score; }
            else { continue; }

            if (new_score < best_score) { best_score = new_score; }

            // accept solution?
            if (new_score < accept_score_) {
                auto msg_copy = std::make_shared<motif_data_structure::MotifStateGraph>(*msg_);
                if (_seen_solution(msg_)) { continue; }
                auto last_n = msg_->last_node()->index();
                auto nj_name = msg_->get_node(nj_)->data()->end_name(ej_);
                auto last_name = msg_->get_node(last_n)->data()->end_name(1);
                msg_copy->add_connection(nj_, last_n, nj_name, last_name);
                return std::make_shared<MotifStateMonteCarloSolutionNew>(msg_copy, new_score);
            }

        }
        std::cout << "stage: " << stage_ << " best_score: " << best_score << " cur_score: " << cur_score << std::endl;

        int j = 0;
        for (int i = org_num_; i < org_num_ + mses_.size(); i++) {
            auto new_ms = mses_[j][rng_.randrange((int) mses_[pos].size() - 1)];
            msg_->replace_state(i, new_ms);
            j++;
        }

        cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);

        step_ = 0;
        stage_ += 1;
    }


    return MotifStateMonteCarloSolutionNewOP(nullptr);
}

void
MotifStateMonteCarlo::start() {
    enumerating_ = true;
}

bool
MotifStateMonteCarlo::finished() {
    if (!enumerating_) {
        throw std::runtime_error("finished() was called but was not enumerating?");
    }

    if (seen_.size() >= max_solutions_) {
        enumerating_ = false;
        return true;
    } else {
        return false;
    }

}


void
MotifStateMonteCarlo::run() {

    auto out = std::ofstream();
    out.open("default.scores");
    out << "design_num,design_score,motif_uses,topology\n";

    auto out_str = std::ofstream();
    out_str.open("default.out");


    auto seen = StringIntMap();

    auto score_num = 0;
    auto pos = 0;
    auto new_ms = motif::MotifStateOP(nullptr);
    auto last_ms = motif::MotifStateOP(nullptr);
    auto cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
    auto new_score = 0.0;
    auto best_score = cur_score;
    auto accept = 0;
    auto motif_used_string = String("");
    mc_.set_temperature(1.0);
    while (stage_ < 100) {
        for (int i = 0; i < 500000; i++) {
            pos = rng_.randrange((int) mses_.size() - 1);
            new_ms = mses_[pos][rng_.randrange((int) mses_[pos].size() - 1)];
            last_ms = msg_->get_node(pos + org_num_)->data()->cur_state;
            msg_->replace_state(pos + org_num_, new_ms);
            new_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
            accept = mc_.accept(cur_score, new_score);
            if (best_score > cur_score) {
                best_score = cur_score;
            }
            if (accept) {
                if (_steric_clash(msg_)) {
                    msg_->replace_state(pos + org_num_, last_ms);
                    continue;
                }

                cur_score = new_score;
            } else {
                msg_->replace_state(pos + org_num_, last_ms);
                continue;
            }
            if (new_score < accept_score_) {
                int j = -1;
                motif_used_string = "";
                for (auto const & n : *msg_) {
                    j++;
                    if (j == 0) { continue; }
                    motif_used_string += n->data()->name() + ";";
                }
                if (seen.find(motif_used_string) != seen.end()) { continue; }
                seen[motif_used_string] = 1;
                if (seen.size() % 10 == 0) {
                    msg_->to_motif_graph()->to_pdb("design." + std::to_string(seen.size() - 1) + ".pdb", 1);
                    std::cout << "solutions: " << seen.size() << std::endl;
                }
                out << score_num << "," << cur_score << "," << motif_used_string << "," << std::endl;
                out_str << msg_->to_motif_graph()->to_str() << std::endl;
            }
        }


        // heat back up
        mc_.set_temperature(10.0);
        for (int i = 0; i < 1000; i++) {
            pos = rng_.randrange((int) mses_.size() - 1);
            new_ms = mses_[pos][rng_.randrange((int) mses_[pos].size() - 1)];
            last_ms = msg_->get_node(pos + org_num_)->data()->cur_state;
            msg_->replace_state(pos + org_num_, new_ms);
            new_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
            accept = mc_.accept(cur_score, new_score);
            if (accept) {
                if (_steric_clash(msg_)) {
                    msg_->replace_state(pos + org_num_, last_ms);
                    continue;
                }
                cur_score = new_score;
            } else {
                msg_->replace_state(pos + org_num_, last_ms);
            }
        }
        mc_.set_temperature(1.0);

        stage_ += 1;
    }
}

double
MotifStateMonteCarlo::get_score(
        structure::BasepairStateOP last_bp) {

    score_ = last_bp->d().distance(end_->d());

    if (!target_an_aligned_end_) {
        r_diff_flip_ = last_bp->r().difference(end_flip_->r());
    } else {
        r_diff_flip_ = last_bp->r().difference(end_->r());
    }
    score_ += r_diff_flip_ * 2;

    return score_;

}

bool
MotifStateMonteCarlo::_steric_clash(
        motif_data_structure::MotifStateGraphOP msg) {

    if (!sterics_) { return false; }

    auto clash = false;
    if (using_lookup_) {
        for (int i = org_num_ + 1; i < msg->size(); i++) {
            for (auto const & b : msg->get_node(i)->data()->cur_state->beads()) {
                clash = lookup_.clash(b);
                if (clash) { return true; }
            }
        }
    }
    return false;
}

float
MotifStateMonteCarlo::perform_motif_swap(
        float cur_score) {
    auto new_ms = motif::MotifStateOP(nullptr);
    auto last_ms = motif::MotifStateOP(nullptr);
    auto pos = rng_.randrange((int) mses_.size() - 1);
    new_ms = mses_[pos][rng_.randrange((int) mses_[pos].size() - 1)];
    last_ms = msg_->get_node(pos + org_num_)->data()->cur_state;
    msg_->replace_state(pos + org_num_, new_ms);
    auto new_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);

    // should we accept?
    auto accept = mc_.accept(cur_score, new_score);

    if (accept) {
        if (_steric_clash(msg_)) {
            msg_->replace_state(pos + org_num_, last_ms);
            return -1;
        }
        return new_score;
    } else {
        msg_->replace_state(pos + org_num_, last_ms);
        return -1;
    }

}


bool
MotifStateMonteCarlo::_seen_solution(
        motif_data_structure::MotifStateGraphOP msg) {

    int j = -1;
    auto motif_used_string = String("");
    for (auto const & n : *msg_) {
        j++;
        if (j == 0) { continue; }
        if (n->data()->name().substr(0, 10) == "HELIX.FLEX") {
            auto spl = base::split_str_by_delimiter(n->data()->name(), ".");
            motif_used_string += spl[0] + "." + spl[1] + "." + spl[2] + ";";
        } else {
            motif_used_string += n->data()->name() + ";";
        }
    }
    if (seen_.find(motif_used_string) != seen_.end()) { return true; }
    else {
        seen_[motif_used_string] = 1;
        return false;
    }
}

}

























