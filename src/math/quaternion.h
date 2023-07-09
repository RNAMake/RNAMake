//
// Created by Joseph Yesselman on 12/16/17.
//

#ifndef RNAMAKE_NEW_QUATERNION_H
#define RNAMAKE_NEW_QUATERNION_H

#define _USE_MATH_DEFINES

#include <math.h>

#include <random>

#include "math/xyz_matrix.h"

namespace math {

class Quaternion {
 public:
  Quaternion() : a_(0), b_(0), c_(0), d_(0) {}

  Quaternion(double val) : a_(val), b_(val), c_(val), d_(val) {}

  Quaternion(double a, double b, double c, double d)
      : a_(a), b_(b), c_(c), d_(d) {}

  Quaternion(Quaternion const &q) : a_(q.a_), b_(q.b_), c_(q.c_), d_(q.d_) {}

 public:
  friend std::ostream &operator<<(std::ostream &, Quaternion const &);

  inline double operator[](int i) const {
    if (i == 0) {
      return a_;
    } else if (i == 1) {
      return b_;
    } else if (i == 2) {
      return c_;
    } else if (i == 3) {
      return d_;
    } else {
      throw std::runtime_error("invalid index");
    }
  }

  inline void operator+=(double const v) {
    a_ += v;
    b_ += v;
    c_ += v;
    d_ += v;
  }

  inline void operator*=(double const v) {
    a_ *= v;
    b_ *= v;
    c_ *= v;
    d_ *= v;
  }

  inline Quaternion operator-() const { return Quaternion(-a_, -b_, -c_, -d_); }

 public:
 public:
  inline double dot(Quaternion const &q) const {
    return a_ * q.a_ + b_ * q.b_ + c_ * q.c_ + d_ * q.d_;
  }

  Matrix get_rotation_matrix();

 public:
  inline double get_a() const { return a_; }

  inline double get_b() const { return b_; }

  inline double get_c() const { return c_; }

  inline double get_d() const { return d_; }

 private:
  double a_, b_, c_, d_;
};

Quaternion get_random_quaternion();

Quaternion get_quaternion_from_matrix(Matrix const &);

inline void dot_vector(std::vector<std::vector<double>> const &m,
                       std::vector<double> const &v, std::vector<double> &vr) {
  for (int i = 0; i < m.size(); i++) {
    vr[i] = 0;
    for (int j = 0; j < m[i].size(); j++) {
      vr[i] += m[i][j] * v[j];
    }
  }
}

inline double norm(std::vector<double> const &v) {
  double norm = 0;
  for (auto const &e : v) {
    norm += e * e;
  }
  return sqrt(norm);
}

void power_iteration(std::vector<std::vector<double>> const &,
                     std::vector<double> &, int);

class AverageQuaternionCalculator {
 public:
  AverageQuaternionCalculator() {
    A_ = std::vector<std::vector<double>>{
        std::vector<double>{0, 0, 0, 0}, std::vector<double>{0, 0, 0, 0},
        std::vector<double>{0, 0, 0, 0}, std::vector<double>{0, 0, 0, 0}};
    qq_ = A_;
    A_scaled_ = A_;

    eigen_values_ = std::vector<double>{0, 0, 0, 0};
    count_ = 0;
  }

  ~AverageQuaternionCalculator() {}

 public:
  void add_quaternion(Quaternion const &q) {
    auto q1 = Quaternion(q);
    if (q1[0] < 0) {
      q1 = -q1;
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        qq_[i][j] = q1[i] * q1[j];
      }
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        A_[i][j] = A_[i][j] + qq_[i][j];
      }
    }
    count_ += 1;
  }

  Quaternion get_average() {
    if (count_ < 1) {
      return Quaternion();
    }

    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        A_scaled_[i][j] = (1 / count_) * A_[i][j];
      }
    }

    power_iteration(A_scaled_, eigen_values_, 100);
    return Quaternion(eigen_values_[0], eigen_values_[1], eigen_values_[2],
                      eigen_values_[3]);
  }

 private:
  double count_;
  std::vector<std::vector<double>> A_scaled_, A_, qq_;
  std::vector<double> eigen_values_;
};

}  // namespace math

#endif  // RNAMAKE_NEW_QUATERNION_H
