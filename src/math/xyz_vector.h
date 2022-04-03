//
//  xyzVector.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__xyzVector__
#define __REDESIGNC__xyzVector__

#include <iostream>
// C++ headers
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// RNAMake Headers
#include <base/string.hpp>

namespace math {

////////////////////////////////////////////////////////////////////////////////
// class declaration                                                          //
////////////////////////////////////////////////////////////////////////////////

class xyzVector {
 private:  // friends //////////////////////////////////////////////////////////
  template <typename>
  friend class xyzMatrix;

 public:  // initiation ////////////////////////////////////////////////////////
  inline xyzVector() = default;

  inline xyzVector(xyzVector const &v) = default;

  inline xyzVector(double const x, double const &y, double const &z)
      : x_(x), y_(y), z_(z) {}

  inline explicit xyzVector(std::vector<double> const &v)
      : x_(v[0]), y_(v[1]), z_(v[2]) {}

 public:  // deletion //////////////////////////////////////////////////////////
  inline ~xyzVector() = default;

 public:  // operators /////////////////////////////////////////////////////////
  /// @brief Copy assignment
  inline xyzVector &operator=(xyzVector const &v) {
    if (this != &v) {
      x_ = v.x_;
      y_ = v.y_;
      z_ = v.z_;
    }
    return *this;
  }
  /// @brief Copy assignment
  template <typename U>
  inline xyzVector &operator=(xyzVector const &v) {
    x_ = v.x_;
    y_ = v.y_;
    z_ = v.z_;
    return *this;
  }
  /// @brief += xyzVector
  template <typename U>
  inline xyzVector &operator+=(xyzVector const &v) {
    x_ += v.x_;
    y_ += v.y_;
    z_ += v.z_;
    return *this;
  }
  /// @brief -= xyzVector
  template <typename U>
  inline xyzVector &operator-=(xyzVector const &v) {
    x_ -= v.x_;
    y_ -= v.y_;
    z_ -= v.z_;
    return *this;
  }
  // @brief = double
  inline xyzVector &operator=(double const &t) {
    x_ = y_ = z_ = t;
    return *this;
  }
  /// @brief += double
  inline xyzVector &operator+=(double const &t) {
    x_ += t;
    y_ += t;
    z_ += t;
    return *this;
  }
  /// @brief -= double
  inline xyzVector &operator-=(double const &t) {
    x_ -= t;
    y_ -= t;
    z_ -= t;
    return *this;
  }
  /// @brief *= double
  inline xyzVector &operator*=(double const &t) {
    x_ *= t;
    y_ *= t;
    z_ *= t;
    return *this;
  }
  /// @brief /= double
  inline xyzVector &operator/=(double const &t) {
    assert(t != double(0));
    double const inv_t(1.0f / t);
    x_ *= inv_t;
    y_ *= inv_t;
    z_ *= inv_t;
    return *this;
  }
  /// @brief -xyzVector (negated copy)
  inline xyzVector operator-() const { return {-x_, -y_, -z_}; }
  /// @brief xyzVector + xyzVector
  friend inline xyzVector operator+(xyzVector const &a, xyzVector const &b) {
    return {a.x_ + b.x_, a.y_ + b.y_, a.z_ + b.z_};
  }
  /// @brief xyzVector + double
  friend inline xyzVector operator+(xyzVector const &v, double const &t) {
    return {v.x_ + t, v.y_ + t, v.z_ + t};
  }
  /// @brief double + xyzVector
  friend inline xyzVector operator+(double const &t, xyzVector const &v) {
    return {t + v.x_, t + v.y_, t + v.z_};
  }
  /// @brief xyzVector - xyzVector
  friend inline xyzVector operator-(xyzVector const &a, xyzVector const &b) {
    return {a.x_ - b.x_, a.y_ - b.y_, a.z_ - b.z_};
  }
  /// @brief xyzVector - double
  friend inline xyzVector operator-(xyzVector const &v, double const &t) {
    return {v.x_ - t, v.y_ - t, v.z_ - t};
  }
  /// @brief double - xyzVector
  friend inline xyzVector operator-(double const &t, xyzVector const &v) {
    return {t - v.x_, t - v.y_, t - v.z_};
  }
  /// @brief xyzVector * double
  friend inline xyzVector operator*(xyzVector const &v, double const &t) {
    return {v.x_ * t, v.y_ * t, v.z_ * t};
  }
  /// @brief double * xyzVector
  friend inline xyzVector operator*(double const &t, xyzVector const &v) {
    return {t * v.x_, t * v.y_, t * v.z_};
  }
  /// @brief xyzVector / double
  friend inline xyzVector operator/(xyzVector const &v, double const &t) {
    assert(t != double(0));
    double const inv_t(double(1) / t);
    return {v.x_ * inv_t, v.y_ * inv_t, v.z_ * inv_t};
  }
  friend inline std::ostream &operator<<(std::ostream &stream,
                                         xyzVector const &v) {
    stream << v.get_str();
    return stream;
  }

 public:  // methods  //////////////////////////////////////////////////////////
  /// @brief Zero
  inline xyzVector &zero() {
    x_ = y_ = z_ = double(0);
    return *this;
  }
  /// @brief Negate
  inline xyzVector &negate() {
    x_ = -x_;
    y_ = -y_;
    z_ = -z_;
    return *this;
  }

  /// @brief Negated copy
  [[nodiscard]] inline xyzVector negated() const { return {-x_, -y_, -z_}; }

  /// @brief Negated: Return via argument (slightly faster)
  inline void negated(xyzVector &a) const {
    a.x_ = -x_;
    a.y_ = -y_;
    a.z_ = -z_;
  }

  /// @brief Normalize
  inline xyzVector &normalize() {
    double const length_ = get_length();
    assert(length_ != double(0));
    double const inv_length(double(1) / length_);
    x_ *= inv_length;
    y_ *= inv_length;
    z_ *= inv_length;
    return *this;
  }

  /// @brief Distance
  [[nodiscard]] inline double distance(xyzVector const &v) const {
    return std::sqrt(square(x_ - v.x_) + square(y_ - v.y_) + square(z_ - v.z_));
  }

  /// @brief Distance squared
  [[nodiscard]] inline double distance_squared(xyzVector const &v) const {
    return square(x_ - v.x_) + square(y_ - v.y_) + square(z_ - v.z_);
  }

  /// @brief Dot product
  [[nodiscard]] inline double dot(xyzVector const &v) const {
    return (x_ * v.x_) + (y_ * v.y_) + (z_ * v.z_);
  }

  /// @brief Cross product
  [[nodiscard]] inline xyzVector cross(xyzVector const &v) const {
    return {(y_ * v.z_) - (z_ * v.y_), (z_ * v.x_) - (x_ * v.z_),
            (x_ * v.y_) - (y_ * v.x_)};
  }

 public:  // Properties: accessors

 public:
  [[nodiscard]] inline String get_str() const {
    return std::to_string(x_) + " " + std::to_string(y_) + " " +
           std::to_string(z_);
  }

  /// @brief double x const
  [[nodiscard]] inline double get_x() const { return x_; }

  /// @brief double y const
  [[nodiscard]] inline double get_y() const { return y_; }

  /// @brief double z const
  [[nodiscard]] inline double get_z() const { return z_; }

  /// @brief Length
  [[nodiscard]] inline double get_length() const {
    return std::sqrt((x_ * x_) + (y_ * y_) + (z_ * z_));
  }

  /// @brief Length squared
  [[nodiscard]] inline double get_length_squared() const {
    return (x_ * x_) + (y_ * y_) + (z_ * z_);
  }

  /// @brief Norm
  [[nodiscard]] inline double get_norm() const {
    return std::sqrt((x_ * x_) + (y_ * y_) + (z_ * z_));
  }

  /// @brief Norm squared
  [[nodiscard]] inline double get_norm_squared() const {
    return (x_ * x_) + (y_ * y_) + (z_ * z_);
  }

  /// @brief Magnitude
  [[nodiscard]] inline double get_magnitude() const {
    return std::sqrt((x_ * x_) + (y_ * y_) + (z_ * z_));
  }

  /// @brief Magnitude squared
  [[nodiscard]] inline double get_magnitude_squared() const {
    return (x_ * x_) + (y_ * y_) + (z_ * z_);
  }

 public:  // Indexers
  /// @brief xyzVector[ i ] const: 0-based index
  inline double const &operator[](int const i) const {
    assert((i >= 0) && (i < 3));
    return (i == 0 ? x_ : (i == 1 ? y_ : z_));
  }

  /// @brief xyzVector[ i ]: 0-based index
  inline double &operator[](int const i) {
    assert((i >= 0) && (i < 3));
    return (i == 0 ? x_ : (i == 1 ? y_ : z_));
  }

  /// @brief xyzVector( i ) const: 1-based index
  inline double const &operator()(int const i) const {
    assert((i > 0) && (i <= 3));
    return (i == 1 ? x_ : (i == 2 ? y_ : z_));
  }

  /// @brief xyzVector( i ): 1-based index
  inline double &operator()(int const i) {
    assert((i > 0) && (i <= 3));
    return (i == 1 ? x_ : (i == 2 ? y_ : z_));
  }

 public:  // Properties: double assignment
  /// @brief x assignment
  inline void set_x(double const &x_a) { x_ = x_a; }

  /// @brief y assignment
  inline void set_y(double const &y_a) { y_ = y_a; }

  /// @brief z assignment
  inline void set_z(double const &z_a) { z_ = z_a; }

 public:  // Comparison
  /// @brief xyzVector == xyzVector
  friend inline bool operator==(xyzVector const &a, xyzVector const &b) {
    return (a.x_ == b.x_) && (a.y_ == b.y_) && (a.z_ == b.z_);
  }

  /// @brief xyzVector != xyzVector
  friend inline bool operator!=(xyzVector const &a, xyzVector const &b) {
    return (a.x_ != b.x_) || (a.y_ != b.y_) || (a.z_ != b.z_);
  }

 private:  // Methods
  /// @brief square( t ) == t * t
  inline static double square(double const &t) { return t * t; }

 private:  // Fields
  /// @brief Coordinates of the 3 coordinate vector
  double x_;
  double y_;
  double z_;

};  // xyzVector

typedef xyzVector Vector;
typedef std::vector<xyzVector> Vectors;

typedef xyzVector Point;
typedef std::vector<Point> Points;

inline Vector vector_from_str(std::string const &s) {
  std::vector<std::string> doubles = base::split_str_by_delimiter(s, " ");
  std::vector<double> point;
  assert(doubles.size() == 3);

  for (auto &i : doubles) {
    point.push_back(std::stod(i.c_str()));
  }

  Vector p(point);
  return p;
}

inline String vector_to_str(Vector const &v) {
  std::stringstream ss;
  ss << v.get_x() << " " << v.get_y() << " " << v.get_z();
  return ss.str();
}

inline String vectors_to_str(Vectors const &vs) {
  std::stringstream ss;
  for (auto const &v : vs) {
    ss << v.get_x() << " " << v.get_y() << " " << v.get_z() << " ";
  }
  return ss.str();
}

inline Vectors vectors_from_str(std::string const &s) {
  std::vector<std::string> doubles = base::split_str_by_delimiter(s, " ");
  std::vector<double> point;
  Vectors vecs;
  for (auto const &d : doubles) {
    point.push_back(std::stod(d.c_str()));
    if (point.size() == 3) {
      Vector vec(point);
      vecs.push_back(vec);
      point = std::vector<double>();
    }
  }

  return vecs;
}

}  // namespace math

#endif /* defined(__REDESIGNC__xyzVector__) */