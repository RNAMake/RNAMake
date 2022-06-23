//
//  BasepairState.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__BasepairState__
#define __REDESIGNC__BasepairState__

#include <iostream>
#include <vector>

// RNAMake Headers
#include "structure/basepair_state.fwd.h"
#include <math/matrix_3x3.hpp>
#include <math/numerical.hpp>
#include <math/transform.hpp>
#include <math/vector_3.hpp>

namespace structure {

class BasepairState {

public:
  inline BasepairState()
      : d_(math::Vector3()), r_(math::Matrix3x3()), r_T_(math::Matrix3x3()),
        sugars_(math::Vector3s(2)), diff_(math::Vector3()),
        diff_sugars_(math::Vector3s(2)) {}

  inline BasepairState(math::Vector3 const &d, math::Matrix3x3 const &r,
                       math::Vector3s const &sugars)
      : d_(d), r_(r), r_T_(math::Matrix3x3(0)), sugars_(sugars),
        diff_(math::Vector3()), diff_sugars_(math::Vector3s(2)) {
    transpose(r_, r_T_);
  }

  inline BasepairState(BasepairState const &b)
      : d_(b.d_), r_(b.r_), r_T_(math::Matrix3x3(0)), sugars_(b.sugars_),
        diff_(math::Vector3()), diff_sugars_(math::Vector3s(2)) {
    transpose(r_, r_T_);
  }

  inline BasepairState(String const &s)
      : r_T_(math::Matrix3x3(0)), diff_(math::Vector3()),
        diff_sugars_(math::Vector3s(2)) {
    auto spl = base::string::split(s, ";");
    if (spl.size() < 3) {
      throw "cannot load BasepairState from String, not the right number of "
            "elements\n";
    }

    d_ = math::Point(spl[0]);
    r_ = math::Matrix(spl[1]);
    sugars_ = math::vectors_from_str(spl[2]);
    transpose(r_, r_T_);
  }

  ~BasepairState() {}

public:
  inline BasepairState copy() const { return BasepairState(d_, r_, sugars_); }

public: // non const methods
  inline void move(math::Vector3 const &p) {
    d_ = d_ + p;
    sugars_[0] = sugars_[0] + p;
    sugars_[1] = sugars_[1] + p;
  }

  inline void transform(math::Matrix3x3 const &r, math::Vector3 const &t,
                        math::Vector3 &dummy) {
    math::dot_vector(r, d_, dummy);
    d_ = dummy + t;

    math::dot_vector(r, sugars_[0], dummy);
    sugars_[0] = dummy + t;

    math::dot_vector(r, sugars_[1], dummy);
    sugars_[1] = dummy + t;

    auto r_new = math::Matrix3x3();
    dot(r_, r, r_new);
    r_ = r_new;
    r_.unitarize();

    transpose(r_, r_T_);
  }

  inline void transform(math::Matrix3x3 const &r, math::Vector3 const &t) {
    auto dummy = math::Vector3();
    transform(r, t, dummy);
  }

public:
  inline void calculate_r_T() { transpose(r_, r_T_); }

  inline void get_transforming_r_and_t(
      BasepairState const
          &o_state, // state with desired rotation and translation
      BasepairState &r_state) {

    // TODO figure out why I have to do this each time???
    calculate_r_T();

    // calculate transforming rotation matrix and store it in r_state (result
    // state)
    dot(r_T_, o_state.r_, r_state.r_);
    r_state.r_.unitarize();
    r_state.calculate_r_T();

    // rotate sugars to new position and store them in r_state
    dot_vectors(r_state.r_T_, o_state.sugars_, r_state.sugars_);

    diff_ = -o_state.d() + d_;
    int i;

    for (i = 0; i < 2; i += 1) {
      r_state.sugars_[i] += diff_;
    }

    for (i = 0; i < 2; i += 1) {
      diff_sugars_[i] = sugars_[i] - r_state.sugars_[i];
    }

    diff_ = (diff_sugars_[0] + diff_sugars_[1]) / 2.0f;
    r_state.d_ = -o_state.d() + diff_ + d_;
  }

  inline void get_transformed_state(BasepairState const &o_state,
                                    BasepairState &r_state) {

    dot(r_, o_state.r_T_, r_state.r_);
    r_state.r_.unitarize();
    math::dot_vector(o_state.r_T_, d_, r_state.d_);
    dot_vectors(o_state.r_T_, sugars_, r_state.sugars_);

    int i;
    for (i = 0; i < 2; i++) {
      r_state.sugars_[i] += o_state.d_;
    }
    r_state.d_ += o_state.d_;
  }

  inline void flip() {
    r_ = transform_1(r_);
    calculate_r_T();
  }

  inline float diff(BasepairStateOP const &state) {
    float diff = d_.distance(state->d());
    diff += _rot_diff(*state) * 2;
    return diff;
  }

  inline float diff(BasepairState const &state) {
    float diff = d_.distance(state.d());
    diff += _rot_diff(state) * 2;
    return diff;
  }

  inline float _rot_diff(BasepairState const &state) {
    float r_diff = r_.difference(state.r());
    // state->flip();
    // float r_diff_2 = r_.difference(state->r());
    // state->flip();
    // if (r_diff > r_diff_2) { r_diff = r_diff_2; }

    return r_diff;
  }

public:
  inline String const to_str() const {
    return d_.get_str() + ";" + r_.get_str() + ";" + sugars_[0].get_str() +
           " " + sugars_[1].get_str();
  }

public: // getters
  inline const math::Vector3 &d() { return d_; }

  inline const math::Vector3 &d() const { return d_; }

  inline const math::Matrix3x3 &r() { return r_; }

  inline const math::Matrix3x3 &r() const { return r_; }

  inline const math::Matrix3x3 &r_T() { return r_T_; }

  inline const math::Matrix3x3 &r_T() const { return r_T_; }

  inline const math::Vector3s &sugars() { return sugars_; }

  inline const math::Vector3s &sugars() const { return sugars_; }

public: // setters
  inline void d(math::Vector3 const &newd) { d_ = newd; }

  inline void r(math::Matrix3x3 const &newr) { r_ = newr; }

  inline void sugars(math::Vector3s const &newsug) { sugars_ = newsug; }

  inline void set(BasepairState const &nstate_) {
    d_ = nstate_.d_;
    r_ = nstate_.r_;
    sugars_ = nstate_.sugars_;
    calculate_r_T();
  }

private:
  math::Vector3 d_;
  math::Matrix3x3 r_;
  math::Matrix3x3 r_T_; // holds the transpose of r_ for speed
  math::Vector3s sugars_;

  /*
  Inclusion of these two tempory variables increases the speed of
  get_transforming_r_and_t() by 240 percent seems worth it for just
  a few more bytes of memory since this a bulk of the calculations
  done in the algorithm -JDY 2014.10.4
  */

  math::Vector3 diff_;         // stores partial products to speed up algorithm
  math::Vector3s diff_sugars_; // stores partial products to speed up algorithm
};

BasepairState str_to_basepairstate(String const &);

BasepairState get_ref_bp_state();

float get_bpstate_rotation_diff(BasepairState const &, BasepairState const &);

std::ostream &operator<<(std::ostream &, const BasepairState &);

inline const float frame_distance(BasepairStateOP const &current,
                                  BasepairStateOP const &end,
                                  BasepairStateOP const &endflip) {

  float score = current->d().distance(end->d());

  float r_diff = end->r().difference(current->r());
  float r_diff_flip = endflip->r().difference(current->r());

  if (r_diff > r_diff_flip) {
    r_diff = r_diff_flip;
  }

  score += r_diff;

  return score;
}

inline const float new_score_function(BasepairStateOP const &current,
                                      BasepairStateOP const &end,
                                      BasepairStateOP const &endflip) {

  float d_diff = current->d().distance(end->d());

  if (d_diff > 25) {
    return d_diff;
  }

  float r_diff = current->r().difference(end->r());
  float r_diff_flip = current->r().difference(endflip->r());

  if (r_diff > r_diff_flip) {
    r_diff = r_diff_flip;
  }

  float scale = (log(150 / d_diff) - 1);
  if (scale > 2) {
    scale = 2;
  }
  if (scale < 0) {
    scale = 0;
  }

  return d_diff + scale * r_diff;
}

inline const float new_score_function_new(BasepairStateOP const &current,
                                          BasepairStateOP const &end,
                                          BasepairStateOP const &endflip) {

  // float d_diff = (current->sugars()[0].distance(end->sugars()[1]) +
  //                 current->sugars()[1].distance(end->sugars()[0]))*0.50;
  float d_diff = current->d().distance(end->d());

  if (d_diff > 25) {
    return d_diff;
  }

  float r_diff = current->r().difference(end->r());
  float r_diff_flip = current->r().difference(endflip->r());

  if (r_diff > r_diff_flip) {
    r_diff = r_diff_flip;
  }

  float scale = (log(150 / d_diff) - 1);
  if (scale > 2) {
    scale = 2;
  }

  return d_diff + scale * r_diff;
}

int are_basepair_states_equal(BasepairState const &, BasepairState const &);

} // namespace structure

#endif /* defined(__REDESIGNC__BasepairState__) */
