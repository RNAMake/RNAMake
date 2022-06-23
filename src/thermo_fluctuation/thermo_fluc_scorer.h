//
//  thermo_fluc_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_scorer__
#define __RNAMake__thermo_fluc_scorer__

#include <stdio.h>

// RNAMake Headers
#include "math/euler.h"
#include "math/hashing.h"
#include "math/matrix_3x3.hpp"
#include "math/vector_3.hpp"
#include "structure/basepair.h"
#include "structure/basepair_state.fwd.h"
#include "structure/basepair_state.h"

namespace thermo_fluctuation {

class ThermoFlucScorer {
public:
  ThermoFlucScorer() {}

  virtual ~ThermoFlucScorer() {}

  virtual inline float score(structure::BasepairStateOP &state_1,
                             structure::BasepairStateOP &state_2) {
    return 0;
  }
};

class FrameScorer : public ThermoFlucScorer {
public:
  FrameScorer()
      : ThermoFlucScorer(), frame_score_(0), r_diff_(0), r_diff_flip_(0) {}

  ~FrameScorer() {}

public:
  inline float score(structure::BasepairStateOP &state_1,
                     structure::BasepairStateOP &state_2) {

    frame_score_ = state_1->d().distance(state_2->d());
    r_diff_ = state_1->r().difference(state_2->r());
    state_2->flip();
    r_diff_flip_ = state_1->r().difference(state_2->r());
    ;
    state_2->flip();

    if (r_diff_ > r_diff_flip_) {
      frame_score_ += r_diff_flip_;
    } else {
      frame_score_ += r_diff_;
    }

    return frame_score_;
  }

private:
  float frame_score_ = 0, r_diff_ = 0, r_diff_flip_ = 0;
};

class FrameScorerDevel : public ThermoFlucScorer {
public:
  FrameScorerDevel()
      : ThermoFlucScorer(), frame_score_(0), r_diff_(0), r_diff_flip_(0),
        weight_d_(1), weight_r_(1) {}

  ~FrameScorerDevel() {}

public:
  inline float score(structure::BasepairStateOP &state_1,
                     structure::BasepairStateOP &state_2) {

    frame_score_ = state_1->d().distance(state_2->d()) * weight_d_;
    r_diff_ = state_1->r().difference(state_2->r());
    state_2->flip();
    r_diff_flip_ = state_1->r().difference(state_2->r());
    state_2->flip();
    if (r_diff_ > r_diff_flip_) {
      frame_score_ += r_diff_flip_ * weight_r_;
    } else {
      frame_score_ += r_diff_ * weight_r_;
    }

    return frame_score_;
  }

public:
  inline void weight_d(float weight_d) { weight_d_ = weight_d; }

  inline void weight_r(float weight_r) { weight_r_ = weight_r; }

public:
  inline float weight_d() { return weight_d_; }

  inline float weight_r() { return weight_r_; }

private:
  float frame_score_ = 0, r_diff_ = 0, r_diff_flip_ = 0;
  float weight_d_ = 0, weight_r_ = 0;
};

class SixDScorer : public ThermoFlucScorer {
public:
  SixDScorer(String const &constraints, structure::BasepairOP ref_bp)
      : ThermoFlucScorer(), ref_bp_(ref_bp) {
    ref_r_t_ = ref_bp_->r().transposed();
    _setup_constraints();
    _parse_constraints(constraints);
  }

public:
  inline float score(structure::BasepairStateOP &state_1,
                     structure::BasepairStateOP &state_2) {

    r1_ = state_1->r();
    d1_ = state_1->d();

    r2_ = state_2->r();
    d2_ = state_2->d();

    dot(ref_r_t_, r1_, rot_);
    rot_t_.unitarize();
    rot_t_ = rot_.transposed();

    dot(r1_, rot_t_, r1_trans_);
    dot(r2_, rot_t_, r2_trans_);

    dot(r1_trans_.transposed(), r2_trans_, r_);
    r_.unitarize();

    auto euler = math::Vector3();
    math::calc_euler(r_, euler);

    d_ = d2_ - d1_;
    for (int i = 0; i < 3; i++) {
      euler[i] = euler[i] * 180 / M_PI;
      euler[i] += 180;
    }

    values_[0] = d_[0];
    values_[1] = d_[1];
    values_[2] = d_[2];
    values_[3] = euler[0];
    values_[4] = euler[1];
    values_[5] = euler[2];

    dist_ = sqrt(values_[0] * values_[0] + values_[1] * values_[1] +
                 values_[2] * values_[2]);
    if (dist_ > 7) {
      return 99;
    }

    fail_ = 0;
    for (int i = 0; i < 6; i++) {
      if (i != 3 && (constraints_[i][0] > values_[i] ||
                     values_[i] > constraints_[i][1])) {
        fail_ = 1;
        break;
      }
      if (i == 3 && (constraints_[i][0] < values_[i] &&
                     values_[i] < constraints_[i][1])) {
        fail_ = 1;
        break;
      }
    }

    if (fail_) {
      return 99;
    } else {
      return 0;
    }
  }

private:
  void _parse_constraints(String const &constraints) {

    auto spl = base::string::split(constraints, ";");
    for (auto const &s : spl) {
      auto spl2 = base::string::split(s, ",");
      if (spl2.size() != 3) {
        throw std::runtime_error("invalid record constraint: " + s);
      }
      auto pos = _parse_constraint_position(spl2[0]);
      auto lower = std::stod(spl2[1]);
      auto upper = std::stod(spl2[2]);
      constraints_[pos] = math::Real2{lower, upper};
    }
  }

  int _parse_constraint_position(String const &constraint_name) {
    if (constraint_name == "x") {
      return 0;
    } else if (constraint_name == "y") {
      return 1;
    } else if (constraint_name == "z") {
      return 2;
    } else if (constraint_name == "a") {
      return 3;
    } else if (constraint_name == "b") {
      return 4;
    } else if (constraint_name == "g") {
      return 5;
    } else {
      throw std::runtime_error("unknown constraint name: " + constraint_name);
    }
  }

  void _setup_constraints() {
    for (int i = 0; i < 3; i++) {
      constraints_[i] = math::Real2{-10, 10};
    }
    for (int i = 3; i < 6; i++) {
      constraints_[i] = math::Real2{0, 360};
    }
  }

private:
  structure::BasepairOP ref_bp_;
  math::Matrix3x3 rot_, rot_t_, ref_r_t_, r_, r1_, r2_, r1_trans_, r2_trans_;
  math::Vector3 d1_, d2_, d_;
  std::array<math::Real2, 6> constraints_;
  math::Real6 values_;
  float dist_;
  bool fail_;
};

typedef std::shared_ptr<ThermoFlucScorer> ThermoFlucScorerOP;

} // namespace thermo_fluctuation

#endif /* defined(__RNAMake__thermo_fluc_scorer__) */
