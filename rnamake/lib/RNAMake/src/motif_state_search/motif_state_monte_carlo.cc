//
// Created by Joseph Yesselman on 11/4/17.
//

#include <fstream>
#include <motif_state_search/motif_state_monte_carlo.h>

void
MotifStateMonteCarlo::setup(
        MotifStateGraphOP msg,
        int ni,
        int nj,
        int ei,
        int ej) {

    msg_ = msg;
    msg_->set_option_value("sterics", false);
    ni_ = ni;
    nj_ = nj;
    ei_ = ei;
    ej_ = ej;

    end_ = msg_->get_node(nj_)->data()->get_end_state(ej_);
    end_flip_ = BasepairStateOP(new BasepairState(end_->copy()));
    end_flip_->flip();

    org_num_ = msg_->size();

    auto i = -1;
    for (auto const & motif_states : mses_) {
        i++;
        auto ms = motif_states[rng_.randrange((int) motif_states.size() - 1)];
        if(i == 0) {
            msg_->add_state(ms, ni_, ei_);
        }
        else {
            msg_->add_state(ms);
        }
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
    auto new_ms = MotifStateOP(nullptr);
    auto last_ms = MotifStateOP(nullptr);
    auto cur_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
    auto new_score = 0.0;
    auto accept = 0;
    auto motif_used_string = String("");
    mc_.set_temperature(1.0);
    for(int stage = 0; stage < 100; stage++ ) {
        std::cout << "stage: " << stage << std::endl;
        for (int i = 0; i < 500000; i++) {
            pos = rng_.randrange((int)mses_.size() - 1);
            new_ms = mses_[pos][rng_.randrange((int) mses_[pos].size() - 1)];
            last_ms = msg_->get_node(pos+org_num_)->data()->cur_state;
            msg_->replace_state(pos+org_num_, new_ms);
            new_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
            accept = mc_.accept(cur_score, new_score);
            if (accept) {
                cur_score = new_score;
            } else {
                msg_->replace_state(pos+org_num_, last_ms);
            }
            if (new_score < 7) {
                int j = -1;
                motif_used_string = "";
                for (auto const & n : *msg_) {
                    j++;
                    if (j == 0) { continue; }
                    motif_used_string += n->data()->name() + ";";
                }
                if(seen.find(motif_used_string) != seen.end()) { continue; }
                seen[motif_used_string] = 1;
                if(seen.size() % 10 == 0) {
                    msg_->to_motif_graph()->to_pdb("design."+std::to_string(seen.size()-1)+".pdb", 1);
                    std::cout << "solutions: " << seen.size() << std::endl;
                }
                out << score_num << "," << cur_score << "," << motif_used_string << "," << std::endl;
                out_str << msg_->to_motif_graph()->to_str() << std::endl;
            }
        }
        // heat back up
        mc_.set_temperature(10.0);
        for(int i = 0; i < 1000; i++)  {
            pos = rng_.randrange((int)mses_.size() - 1);
            new_ms = mses_[pos][rng_.randrange((int) mses_[pos].size() - 1)];
            last_ms = msg_->get_node(pos+org_num_)->data()->cur_state;
            msg_->replace_state(pos+org_num_, new_ms);
            new_score = get_score(msg_->last_node()->data()->cur_state->end_states()[1]);
            accept = mc_.accept(cur_score, new_score);
            if (accept) {
                cur_score = new_score;
            } else {
                msg_->replace_state(pos+org_num_, last_ms);
            }
        }
        mc_.set_temperature(1.0);
    }
}

double
MotifStateMonteCarlo::get_score(
        BasepairStateOP last_bp) {

    score_ = last_bp->d().distance(end_->d());

    //r_diff_       = last_bp->r().difference(end_->r());
    r_diff_flip_  = last_bp->r().difference(end_flip_->r());
    score_ += r_diff_flip_*2;

    //if( r_diff_ > r_diff_flip_) { score_ += 2*r_diff_flip_; }
    //else                        { score_ += 2*r_diff_;      }

    return score_;

}
