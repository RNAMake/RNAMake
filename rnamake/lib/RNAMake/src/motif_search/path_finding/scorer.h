//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SCORER_H
#define RNAMAKE_NEW_PATH_FINDING_SCORER_H

#include <structure/basepair_state.h>
#include <motif_search/path_finding/node.h>

namespace motif_search {
namespace path_finding {

class Scorer {
public:
    Scorer() {}

    virtual
    ~Scorer() {}

    virtual
    Scorer *
    clone() const = 0;

public:
    void
    set_target(
            structure::BasepairStateOP target,
            bool target_an_aligned_end) {
        target_ = target;
        target_an_aligned_end_ = target_an_aligned_end;
        target_flip_ = std::make_shared<structure::BasepairState>(target_->copy());
        target_flip_->flip();
    }

public:
    virtual
    inline
    float
    score(
            Node const &) = 0;

    virtual
    inline
    float
    score(
            motif::MotifState &,
            Node const &) = 0;

    virtual
    inline
    float
    accept_score(Node const & node)  {
        best_score_ = 1000;
        int i = -1;
        for (auto const & state : node.state()->end_states()) {
            i++;
            if (i == 0) { continue; }

            score_ = state->d().distance(target_->d());

            if (target_an_aligned_end_) { r_diff_ = state->r().difference(target_->r());      }
            else                        { r_diff_ = state->r().difference(target_flip_->r()); }
            score_ += 2 * r_diff_;

            if (score_ < best_score_) { best_score_ = score_; }
        }
        return best_score_;
    }

    virtual
    void
    set_dummy(float dummy) {}


protected:

    inline
    float
    _weighted_score(
            structure::BasepairStateOP const & current,
            structure::BasepairStateOP const & end,
            structure::BasepairStateOP const & endflip) {
        d_diff_= current->d().distance(end->d());
        if (d_diff_ > 25) { return d_diff_; }

        if (target_an_aligned_end_) {  r_diff_ = current->r().difference(target_->r());      }
        else                        {  r_diff_ = current->r().difference(target_flip_->r()); }

        scale_ = (log(150 / d_diff_) - 1);
        if (scale_ > 2) { scale_ = 2; }

        return d_diff_ + scale_*r_diff_;
    }


protected:
    structure::BasepairStateOP target_, target_flip_;
    float best_score_, score_, r_diff_, r_diff_flip_, d_diff_, scale_;
    bool target_an_aligned_end_;
};

typedef std::shared_ptr<Scorer> ScorerOP;

class GreedyBestFirst : public Scorer {
public:
    GreedyBestFirst() : Scorer() {}

    Scorer *
    clone() const { return new GreedyBestFirst(*this); };

public:
    inline
    float
    score(
            Node const & node) {
        best_score_ = 1000;int i = -1;
        for (auto const & state : node.state()->end_states()) {
            i++;
            if (i == 0) { continue; }

            score_ = _weighted_score(state, target_, target_flip_);

            if (score_ < best_score_) {
                best_score_ = score_;
            }
        }
        return best_score_;
    }

    inline
    float
    score(
            motif::MotifState & ms,
            Node const & node) {
        best_score_ = 1000;
        int i = -1;
        for (auto const & state : ms.end_states()) {
            i++;
            if (i == 0) { continue; }

            score_ = _weighted_score(state, target_, target_flip_);

            if (score_ < best_score_) {
                best_score_ = score_;
            }

        }
        return best_score_;

    }

};


}
}

#endif //RNAMAKE_NEW_PATH_FINDING_SCORER_H