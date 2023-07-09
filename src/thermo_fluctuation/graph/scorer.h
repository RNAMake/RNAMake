//
// Created by Joseph Yesselman on 2019-04-09.
//

#ifndef RNAMAKE_NEW_THERMO_FLUC_GRAPH_SCORER_H
#define RNAMAKE_NEW_THERMO_FLUC_GRAPH_SCORER_H

#include <structure/basepair_state.h>

namespace thermo_fluctuation {
namespace graph {

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
    setup(
            bool target_an_aligned_end) {
        target_an_aligned_end_ = target_an_aligned_end;
    }

public:
    virtual
    inline
    float
    score(
            structure::BasepairState const &,
            structure::BasepairState const &) = 0;

protected:
    bool target_an_aligned_end_;
};

typedef std::shared_ptr<Scorer> ScorerOP;

class FrameScorer : public Scorer {
public:
    FrameScorer() : Scorer() {}

    ~FrameScorer() {}

    Scorer *
    clone() const { return new FrameScorer(*this); };

public:
    inline
    float
    score(
            structure::BasepairState const & state_1,
            structure::BasepairState const & state_2) {

        score_ = state_1.d().distance(state_2.d());


        if(target_an_aligned_end_) {
            r_diff_ = state_1.r().difference(state_2.r());
        }
        else {
            flipped_ = state_2.r().get_flip_orientation();
            r_diff_ = state_1.r().difference(flipped_);
        }
        //return score_ + 2*r_diff_;
        return score_ + 1*r_diff_;

    }

private:
    float r_diff_, score_;
    math::Matrix flipped_;

};

// copy old code scoring
class OldFrameScorer : public Scorer {
public:
    OldFrameScorer() : thermo_fluctuation::graph::Scorer() {}

    ~OldFrameScorer() {}

    thermo_fluctuation::graph::Scorer *
    clone() const { return new OldFrameScorer(*this); };

public:

    inline
    float
    score(
            structure::BasepairState const & state_1,
            structure::BasepairState const & state_2) {
        score_ = state_1.d().distance(state_2.d());
        r_diff_ = state_1.r().difference(state_2.r());
        flipped_ = state_2.r().get_flip_orientation();
        r_diff_flip_ = state_1.r().difference(flipped_);
        if (r_diff_ > r_diff_flip_) { score_ += r_diff_flip_; }
        else { score_ += r_diff_; }
        return score_;
    }

private:
    float r_diff_, r_diff_flip_,  score_;
    math::Matrix flipped_;
};

}
}


#endif //RNAMAKE_NEW_THERMO_FLUC_GRAPH_SCORER_H
