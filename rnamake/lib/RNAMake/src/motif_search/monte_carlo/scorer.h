//
// Created by Joseph Yesselman on 3/23/19.
//

#ifndef RNAMAKE_NEW_MONTE_CARLO_SCORER_H
#define RNAMAKE_NEW_MONTE_CARLO_SCORER_H

#include <structure/basepair_state.h>
#include <motif/motif_state.h>

namespace motif_search {
namespace monte_carlo {

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
            structure::BasepairState const & bps) = 0;

    virtual
    inline
    float
    accept_score(
            structure::BasepairState const & bps)  {
        int i = -1;

        score_ = bps.d().distance(target_->d());

        if (target_an_aligned_end_) { r_diff_ = bps.r().difference(target_->r());      }
        else                        { r_diff_ = bps.r().difference(target_flip_->r()); }
        score_ += 2 * r_diff_;

        return score_;
    }


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

class GreedyScorer : public Scorer {
public:
    GreedyScorer() : Scorer() {}

    Scorer *
    clone() const { return new GreedyScorer(*this); };

public:
    inline
    float
    score(
            structure::BasepairState const & bps) {
        return accept_score(bps);
    }

};


}
}

#endif //RNAMAKE_NEW_MONTE_CARLO_SCORER_H
