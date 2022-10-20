//
// Created by Joseph Yesselman on 3/17/19.
//

#ifndef RNAMAKE_NEW_PATH_FINDING_SCORER_H
#define RNAMAKE_NEW_PATH_FINDING_SCORER_H

#include <motif_search/path_finding/node.h>
#include <structure/basepair_state.h>

namespace motif_search {
namespace path_finding {

class Scorer {
public:
  Scorer() {}

  virtual ~Scorer() {}

  virtual Scorer *clone() const = 0;

public:
  void set_target(structure::BasepairStateOP target,
                  bool target_an_aligned_end) {
    target_ = target;
    target_an_aligned_end_ = target_an_aligned_end;
    target_flip_ = std::make_shared<structure::BasepairState>(target_->copy());
    target_flip_->flip();
  }

public:
  virtual inline float score(Node const &) = 0;

  virtual inline float score(motif::MotifState &, Node const &) = 0;

  virtual inline float accept_score(Node const &node) {
    best_score_ = 1000;
    int i = -1;
    for (auto const &state : node.state()->end_states()) {
      i++;
      if (i == 0) {
        continue;
      }

      score_ = state->d().distance(target_->d());

      if (target_an_aligned_end_) {
        r_diff_ = state->r().difference(target_->r());
      } else {
        r_diff_ = state->r().difference(target_flip_->r());
      }
      score_ += 2 * r_diff_;

      if (score_ < best_score_) {
        best_score_ = score_;
      }
    }
    return best_score_;
  }

  virtual void set_dummy(float dummy) {}

protected:
  inline float _weighted_score(structure::BasepairStateOP const &current,
                               structure::BasepairStateOP const &end,
                               structure::BasepairStateOP const &endflip) {
    d_diff_ = current->d().distance(end->d());
    if (d_diff_ > 25) {
      return d_diff_;
    }

    if (target_an_aligned_end_) {
      r_diff_ = current->r().difference(target_->r());
    } else {
      r_diff_ = current->r().difference(target_flip_->r());
    }

    scale_ = (log(150 / d_diff_) - 1);
    if (scale_ > 2) {
      scale_ = 2;
    }

    return d_diff_ + scale_ * r_diff_;
  }

protected:
  structure::BasepairStateOP _target, _target_flip;
  float _best_score, _score, _r_diff, _r_diff_flip, _d_diff, _scale;
  bool _target_an_aligned_end;
};

typedef std::shared_ptr<Scorer> ScorerOP;

class GreedyScorer : public Scorer {
public:
  GreedyScorer() : Scorer() {}

  Scorer *clone() const { return new GreedyScorer(*this); };

public:
  inline float score(Node const &node) {
    _best_score = 1000;
    int i = -1;
    for (auto const &state : node.state()->end_states()) {
      i++;
      if (i == 0) {
        continue;
      }

      _score = _weighted_score(state, _target, _target_flip);

      if (_score < _best_score) {
        _best_score = _score;
      }
    }
    return _best_score;
  }

  inline float score(motif::MotifState &ms, Node const &node) {
    _best_score = 1000;
    int i = -1;
    for (auto const &state : ms.end_states()) {
      i++;
      if (i == 0) {
        continue;
      }

      _score = _weighted_score(state, _target, _target_flip);

      if (_score < _best_score) {
        _best_score = _score;
      }
    }
    return _best_score;
  }
};

class AstarScorer : public Scorer {
public:
  AstarScorer() : Scorer() {
    _g = 0;
    _g = 0;
    _ss_score_weight = 0.10;
    _level_weight = 3;
  }

  Scorer *clone() const { return new AstarScorer(*this); };

public:
  inline float score(Node const &node) {
    _best_score = 1000;
    int i = -1;
    for (auto const &state : node.state()->end_states()) {
      i++;
      if (i == 0) {
        continue;
      }

      _score = _weighted_score(state, _target, _target_flip);

      if (_score < _best_score) {
        _best_score = _score;
      }
    }
    _h = _best_score;
    _g = node.ss_score() * _ss_score_weight;
    if (node.level() > 3) {
      _g += node.level() * _level_weight;
    }
    std::cout << _h << " " << _g << std::endl;
    return _g + _h;
  }

  inline float score(motif::MotifState &ms, Node const &node) {
    _best_score = 1000;
    int i = -1;
    for (auto const &state : ms.end_states()) {
      i++;
      if (i == 0) {
        continue;
      }

      _score = _weighted_score(state, _target, _target_flip);

      if (_score < _best_score) {
        _best_score = _score;
      }
    }

    _h = _best_score;
    _g = (node.ss_score() + ms.score()) * _ss_score_weight;
    if (node.level() + 1 > 3) {
      _g += (node.level() + 1) * _level_weight;
    }
    return _h + _g;
  }

private:
  float _g, _h;
  float _ss_score_weight, _level_weight;
};

} // namespace path_finding
} // namespace motif_search

#endif // RNAMAKE_NEW_PATH_FINDING_SCORER_H
