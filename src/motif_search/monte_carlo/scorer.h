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


protected:
    structure::BasepairStateOP target_, target_flip_;
    float best_score_, score_, r_diff_, r_diff_flip_, d_diff_, scale_;
    bool target_an_aligned_end_;
};

typedef std::shared_ptr<Scorer> ScorerOP;

class DefaultScorer : public Scorer {
public:
    DefaultScorer() : Scorer() {}

    Scorer *
    clone() const { return new DefaultScorer(*this); };

public:
    inline
    float
    score(
            structure::BasepairState const & bps) {
        score_ = bps.d().distance(target_->d());

        if (target_an_aligned_end_) { r_diff_ = bps.r().difference(target_->r());      }
        else                        { r_diff_ = bps.r().difference(target_flip_->r()); }
        score_ += 5 * r_diff_;
        //score_ += r_diff_;
        //score_ = bps.sugars()[0].distance(target_->sugars()[1]) + bps.sugars()[1].distance(target_->sugars()[0]);
        return score_;
    }
};

class ScaledScorer : public Scorer {
public:
    ScaledScorer(
            float d_scale,
            float r_scale) : Scorer() {}

    Scorer *
    clone() const { return new ScaledScorer(*this); };

public:
    inline
    float
    score(
            structure::BasepairState const & bps) {
        score_ = bps.d().distance(target_->d())*d_scale_;

        if (target_an_aligned_end_) { r_diff_ = bps.r().difference(target_->r());      }
        else                        { r_diff_ = bps.r().difference(target_flip_->r()); }
        score_ += r_scale_ * r_diff_;
        //score_ += r_diff_;
        //score_ = bps.sugars()[0].distance(target_->sugars()[1]) + bps.sugars()[1].distance(target_->sugars()[0]);
        return score_;
    }

private:
    float d_scale_, r_scale_;


};

class ScorerFactory  {
public:
    ScorerFactory() {}

public:
    ScorerOP
    get_scorer(
            String const & scorer_name) {
        if(scorer_name == "default") {
            return ScorerOP(std::make_shared<DefaultScorer>());
        }
        else if(scorer_name == "sugar_dist") {
            return ScorerOP(std::make_shared<DefaultScorer>());
        }
        else {
            throw std::runtime_error(scorer_name + " is not a valid scorer name");
        }

    }

};

}
}

#endif //RNAMAKE_NEW_MONTE_CARLO_SCORER_H
