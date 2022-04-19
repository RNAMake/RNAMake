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
// std headers
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
      : _x(x), _y(y), _z(z) {}

  inline explicit xyzVector(std::vector<double> const &v)
      : _x(v[0]), _y(v[1]), _z(v[2]) {}

 public:  // deletion //////////////////////////////////////////////////////////
  inline ~xyzVector() = default;

 public:  // operators /////////////////////////////////////////////////////////
  /// @brief Copy assignment
  inline xyzVector &operator=(xyzVector const &v) {
    if (this != &v) {
      _x = v._x;
      _y = v._y;
      _z = v._z;
    }
    return *this;
  }
  /// @brief += xyzVector
  inline xyzVector &operator+=(xyzVector const &v) {
    _x += v._x;
    _y += v._y;
    _z += v._z;
    return *this;
  }
  /// @brief -= xyzVector
  inline xyzVector &operator-=(xyzVector const &v) {
    _x -= v._x;
    _y -= v._y;
    _z -= v._z;
    return *this;
  }
  // @brief = double
  inline xyzVector &operator=(double const &t) {
    _x = _y = _z = t;
    return *this;
  }
  /// @brief += double
  inline xyzVector &operator+=(double const &t) {
    _x += t;
    _y += t;
    _z += t;
    return *this;
  }
  /// @brief -= double
  inline xyzVector &operator-=(double const &t) {
    _x -= t;
    _y -= t;
    _z -= t;
    return *this;
  }
  /// @brief *= double
  inline xyzVector &operator*=(double const &t) {
    _x *= t;
    _y *= t;
    _z *= t;
    return *this;
  }
  /// @brief /= double
  inline xyzVector &operator/=(double const &t) {
    assert(t != double(0));
    double const inv_t(1.0f / t);
    _x *= inv_t;
    _y *= inv_t;
    _z *= inv_t;
    return *this;
  }
  /// @brief -xyzVector (negated copy)
  inline xyzVector operator-() const { return {-_x, -_y, -_z}; }
  /// @brief xyzVector + xyzVector
  friend inline xyzVector operator+(xyzVector const &a, xyzVector const &b) {
    return {a._x + b._x, a._y + b._y, a._z + b._z};
  }
  /// @brief xyzVector + double
  friend inline xyzVector operator+(xyzVector const &v, double const &t) {
    return {v._x + t, v._y + t, v._z + t};
  }
  /// @brief double + xyzVector
  friend inline xyzVector operator+(double const &t, xyzVector const &v) {
    return {t + v._x, t + v._y, t + v._z};
  }
  /// @brief xyzVector - xyzVector
  friend inline xyzVector operator-(xyzVector const &a, xyzVector const &b) {
    return {a._x - b._x, a._y - b._y, a._z - b._z};
  }
  /// @brief xyzVector - double
  friend inline xyzVector operator-(xyzVector const &v, double const &t) {
    return {v._x - t, v._y - t, v._z - t};
  }
  /// @brief double - xyzVector
  friend inline xyzVector operator-(double const &t, xyzVector const &v) {
    return {t - v._x, t - v._y, t - v._z};
  }
  /// @brief xyzVector * double
  friend inline xyzVector operator*(xyzVector const &v, double const &t) {
    return {v._x * t, v._y * t, v._z * t};
  }
  /// @brief double * xyzVector
  friend inline xyzVector operator*(double const &t, xyzVector const &v) {
    return {t * v._x, t * v._y, t * v._z};
  }
  /// @brief xyzVector / double
  friend inline xyzVector operator/(xyzVector const &v, double const &t) {
    assert(t != double(0));
    double const inv_t(double(1) / t);
    return {v._x * inv_t, v._y * inv_t, v._z * inv_t};
  }
  friend inline std::ostream &operator<<(std::ostream &stream,
                                         xyzVector const &v) {
    stream << v.get_str();
    return stream;
  }

 public:  // methods  //////////////////////////////////////////////////////////
  /// @brief Zero
  inline xyzVector &zero() {
    _x = _y = _z = double(0);
    return *this;
  }
  /// @brief Negate
  inline xyzVector &negate() {
    _x = -_x;
    _y = -_y;
    _z = -_z;
    return *this;
  }

  /// @brief Negated copy
  [[nodiscard]] inline xyzVector negated() const { return {-_x, -_y, -_z}; }

  /// @brief Negated: Return via argument (slightly faster)
  inline void negated(xyzVector &a) const {
    a._x = -_x;
    a._y = -_y;
    a._z = -_z;
  }

  /// @brief Normalize
  inline xyzVector &normalize() {
    double const length_ = get_length();
    assert(length_ != double(0));
    double const inv_length(double(1) / length_);
    _x *= inv_length;
    _y *= inv_length;
    _z *= inv_length;
    return *this;
  }

  /// @brief Distance
  [[nodiscard]] inline double distance(xyzVector const &v) const {
    return std::sqrt(square(_x - v._x) + square(_y - v._y) + square(_z - v._z));
  }

  /// @brief Distance squared
  [[nodiscard]] inline double distance_squared(xyzVector const &v) const {
    return square(_x - v._x) + square(_y - v._y) + square(_z - v._z);
  }

  /// @brief Dot product
  [[nodiscard]] inline double dot(xyzVector const &v) const {
    return (_x * v._x) + (_y * v._y) + (_z * v._z);
  }

  /// @brief Cross product
  [[nodiscard]] inline xyzVector cross(xyzVector const &v) const {
    return {(_y * v._z) - (_z * v._y), (_z * v._x) - (_x * v._z),
            (_x * v._y) - (_y * v._x)};
  }

 public:  // Properties: accessors

 public:
  [[nodiscard]] inline String get_str() const {
    return std::to_string(_x) + " " + std::to_string(_y) + " " +
           std::to_string(_z);
  }

  /// @brief get x value
  [[nodiscard]] inline double get_x() const { return _x; }

  /// @brief get y value
  [[nodiscard]] inline double get_y() const { return _y; }

  /// @brief get z value
  [[nodiscard]] inline double get_z() const { return _z; }

  /// @brief Length
  [[nodiscard]] inline double get_length() const {
    return std::sqrt((_x * _x) + (_y * _y) + (_z * _z));
  }

  /// @brief Length squared
  [[nodiscard]] inline double get_length_squared() const {
    return (_x * _x) + (_y * _y) + (_z * _z);
  }

  /// @brief Norm
  [[nodiscard]] inline double get_norm() const {
    return std::sqrt((_x * _x) + (_y * _y) + (_z * _z));
  }

  /// @brief Norm squared
  [[nodiscard]] inline double get_norm_squared() const {
    return (_x * _x) + (_y * _y) + (_z * _z);
  }

  /// @brief Magnitude
  [[nodiscard]] inline double get_magnitude() const {
    return std::sqrt((_x * _x) + (_y * _y) + (_z * _z));
  }

  /// @brief Magnitude squared
  [[nodiscard]] inline double get_magnitude_squared() const {
    return (_x * _x) + (_y * _y) + (_z * _z);
  }

 public:  // Indexers
  /// @brief xyzVector[ i ] const: 0-based index
  inline double const &operator[](int const i) const {
    assert((i >= 0) && (i < 3));
    return (i == 0 ? _x : (i == 1 ? _y : _z));
  }

  /// @brief xyzVector[ i ]: 0-based index
  inline double &operator[](int const i) {
    assert((i >= 0) && (i < 3));
    return (i == 0 ? _x : (i == 1 ? _y : _z));
  }
  
 public:  // Properties: double assignment
  /// @brief x assignment
  inline void set_x(double const &x_a) { _x = x_a; }

  /// @brief y assignment
  inline void set_y(double const &y_a) { _y = y_a; }

  /// @brief z assignment
  inline void set_z(double const &z_a) { _z = z_a; }

 public:  // Comparison
  /// @brief xyzVector == xyzVector
  friend inline bool operator==(xyzVector const &a, xyzVector const &b) {
    return (a._x == b._x) && (a._y == b._y) && (a._z == b._z);
  }

  /// @brief xyzVector != xyzVector
  friend inline bool operator!=(xyzVector const &a, xyzVector const &b) {
    return (a._x != b._x) || (a._y != b._y) || (a._z != b._z);
  }

 private:  // Methods
  /// @brief square( t ) == t * t
  inline static double square(double const &t) { return t * t; }

 private:  // Fields
  /// @brief Coordinates of the 3 coordinate vector
  double _x;
  double _y;
  double _z;

};  // xyzVector

typedef xyzVector Vector;
typedef std::vector<xyzVector> Vectors;

typedef xyzVector Point;
typedef std::vector<Point> Points;

inline Vector vector_from_str(String const & s) {
  Strings doubles = base::string::split(s, " ");
  Reals point;
  assert(doubles.size() == 3);

  for (auto &i : doubles) {
    point.push_back(std::stod(i.c_str()));
  }

  Vector p(point);
  return p;
}

inline Vectors vectors_from_str(String const &s) {
  Strings doubles = base::string::split(s, " ");
  Reals point;
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