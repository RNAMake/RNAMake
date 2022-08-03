//
// Created by Joseph Yesselman on 2019-03-30.
//

#ifndef RNAMAKE_NEW_EXHAUSTIVE_SCORER_H
#define RNAMAKE_NEW_EXHAUSTIVE_SCORER_H

#include <motif/motif_state.h>
#include <structure/basepair_state.h>
#include "base/types.hpp"

namespace motif_search {
namespace exhaustive {

class Scorer {
public:
  Scorer() {}

  virtual ~Scorer() {}

  virtual Scorer *clone() const = 0;

public:
  void set_target(structure::BasepairStateOP target,
                  bool target_an_aligned_end) {
    _target = target;
    _target_an_aligned_end = target_an_aligned_end;
    _target_flip = std::make_shared<structure::BasepairState>(_target->copy());
    _target_flip->flip();
  }

public:
  virtual inline float score(structure::BasepairState const &bps) = 0;

protected:
  structure::BasepairStateOP _target, _target_flip;
  float _best_score, _score, _r_diff, _r_diff_flip, _d_diff, _scale;
  bool _target_an_aligned_end;
};

typedef std::shared_ptr<Scorer> ScorerOP;

class DefaultScorer : public Scorer {
public:
  DefaultScorer() : Scorer() {}

  Scorer *clone() const { return new DefaultScorer(*this); };

public:
  inline float score(structure::BasepairState const &bps) {
    _score = bps.d().distance(_target->d());

    if (_target_an_aligned_end) {
      _r_diff = bps.r().difference(_target->r());
    } else {
      _r_diff = bps.r().difference(_target_flip->r());
    }
    _score += 2 * _r_diff;
    // score_ += r_diff_;
    // score_ = bps.sugars()[0].distance(target_->sugars()[1]) +
    // bps.sugars()[1].distance(target_->sugars()[0]);
    return _score;
  }
};

class ScorerFactory {
public:
  ScorerFactory() {}

public:
  ScorerOP get_scorer(String const &scorer_name) {
    if (scorer_name == "default") {
      return ScorerOP(std::make_shared<DefaultScorer>());
    } else if (scorer_name == "sugar_dist") {
      return ScorerOP(std::make_shared<DefaultScorer>());
    } else {
      throw std::runtime_error(scorer_name + " is not a valid scorer name");
    }
  }
};

} // namespace exhaustive

} // namespace motif_search

#endif // RNAMAKE_NEW_EXHAUSTIVE_SCORER_H
