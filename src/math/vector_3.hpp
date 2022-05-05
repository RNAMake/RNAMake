//
//  Vector3.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Vector3__
#define __REDESIGNC__Vector3__

#include <iostream>
// std headers
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// RNAMake Headers
#include <base/exception.hpp>
#include <base/string.hpp>

namespace math {

////////////////////////////////////////////////////////////////////////////////
// class declaration                                                          //
////////////////////////////////////////////////////////////////////////////////

class Vector3 {
private: // friends //////////////////////////////////////////////////////////
  template <typename> friend class xyzMatrix;

public: // initiation ////////////////////////////////////////////////////////
  inline Vector3() = default;
  inline Vector3(const Vector3 &v) = default;
  inline Vector3(const double x, const double &y, const double &z)
      : _x(x), _y(y), _z(z) {}
  inline explicit Vector3(const Reals &v) : _x(v[0]), _y(v[1]), _z(v[2]) {}

public: // deletion //////////////////////////////////////////////////////////
  inline ~Vector3() = default;

public: // operators /////////////////////////////////////////////////////////
  /// @brief Copy assignment
  inline Vector3 &operator=(const Vector3 &v) {
    if (this != &v) {
      _x = v._x;
      _y = v._y;
      _z = v._z;
    }
    return *this;
  }
  // @brief = double
  inline Vector3 &operator=(const double &t) {
    _x = _y = _z = t;
    return *this;
  }
  /// @brief -Vector3 (negated copy)
  inline Vector3 operator-() const { return {-_x, -_y, -_z}; }
  /// @brief Vector3 + Vector3
  friend inline Vector3 operator+(const Vector3 &a, const Vector3 &b) {
    return {a._x + b._x, a._y + b._y, a._z + b._z};
  }
  /// @brief Vector3 + double
  friend inline Vector3 operator+(const Vector3 &v, const double &t) {
    return {v._x + t, v._y + t, v._z + t};
  }
  /// @brief double + Vector3
  friend inline Vector3 operator+(const double &t, const Vector3 &v) {
    return {t + v._x, t + v._y, t + v._z};
  }
  /// @brief Vector3 - Vector3
  friend inline Vector3 operator-(const Vector3 &a, const Vector3 &b) {
    return {a._x - b._x, a._y - b._y, a._z - b._z};
  }
  /// @brief Vector3 - double
  friend inline Vector3 operator-(const Vector3 &v, const double &t) {
    return {v._x - t, v._y - t, v._z - t};
  }
  /// @brief double - Vector3
  friend inline Vector3 operator-(const double &t, const Vector3 &v) {
    return {t - v._x, t - v._y, t - v._z};
  }
  /// @brief Vector3 * double
  friend inline Vector3 operator*(const Vector3 &v, const double &t) {
    return {v._x * t, v._y * t, v._z * t};
  }
  /// @brief double * Vector3
  friend inline Vector3 operator*(const double &t, const Vector3 &v) {
    return {t * v._x, t * v._y, t * v._z};
  }
  /// @brief Vector3 / double
  friend inline Vector3 operator/(const Vector3 &v, const double &t) {
    if (t == 0) {
      String msg = "trying to divide a vector by zero!";
      base::log_and_throw<base::MathException>(msg);
    }
    double const inv_t(double(1) / t);
    return {v._x * inv_t, v._y * inv_t, v._z * inv_t};
  }
  /// @brief Vector3 <<
  friend inline std::ostream &operator<<(std::ostream &stream,
                                         const Vector3 &v) {
    stream << v.get_str();
    return stream;
  }
  /// @brief += Vector3
  inline Vector3 &operator+=(const Vector3 &v) {
    _x += v._x;
    _y += v._y;
    _z += v._z;
    return *this;
  }
  /// @brief += double
  inline Vector3 &operator+=(const double &t) {
    _x += t;
    _y += t;
    _z += t;
    return *this;
  }
  /// @brief -= Vector3
  inline Vector3 &operator-=(const Vector3 &v) {
    _x -= v._x;
    _y -= v._y;
    _z -= v._z;
    return *this;
  }
  /// @brief -= double
  inline Vector3 &operator-=(const double &t) {
    _x -= t;
    _y -= t;
    _z -= t;
    return *this;
  }
  /// @brief *= double
  inline Vector3 &operator*=(const double &t) {
    _x *= t;
    _y *= t;
    _z *= t;
    return *this;
  }
  /// @brief /= double
  inline Vector3 &operator/=(const double &t) {
    if (t == 0) {
      String msg = "trying to divide a vector by zero!";
      base::log_and_throw<base::MathException>(msg);
    }
    const double inv_t(1.0f / t);
    _x *= inv_t;
    _y *= inv_t;
    _z *= inv_t;
    return *this;
  }

public: // comparison ///////////////////////////////////////////////////////
  /// @brief Vector3 == Vector3
  friend inline bool operator==(const Vector3 &a, const Vector3 &b) {
    return (a._x == b._x) && (a._y == b._y) && (a._z == b._z);
  }
  /// @brief Vector3 != Vector3
  friend inline bool operator!=(const Vector3 &a, const Vector3 &b) {
    return (a._x != b._x) || (a._y != b._y) || (a._z != b._z);
  }

public: // methods  //////////////////////////////////////////////////////////
  /// @brief Cross product
  [[nodiscard]] inline Vector3 cross(const Vector3 &v) const {
    return {(_y * v._z) - (_z * v._y), (_z * v._x) - (_x * v._z),
            (_x * v._y) - (_y * v._x)};
  }
  /// @brief Distance
  [[nodiscard]] inline double distance(const Vector3 &v) const {
    return std::sqrt(square(_x - v._x) + square(_y - v._y) + square(_z - v._z));
  }
  /// @brief Distance squared
  [[nodiscard]] inline double distance_squared(const Vector3 &v) const {
    return square(_x - v._x) + square(_y - v._y) + square(_z - v._z);
  }
  /// @brief Dot product
  [[nodiscard]] inline double dot(const Vector3 &v) const {
    return (_x * v._x) + (_y * v._y) + (_z * v._z);
  }

  /// @brief Zero
  inline Vector3 &zero() {
    _x = _y = _z = double(0);
    return *this;
  }
  /// @brief Negate
  inline Vector3 &negate() {
    _x = -_x;
    _y = -_y;
    _z = -_z;
    return *this;
  }

  /// @brief Negated copy
  [[nodiscard]] inline Vector3 negated() const { return {-_x, -_y, -_z}; }

  /// @brief Negated: Return via argument (slightly faster)
  inline void negated(Vector3 &a) const {
    a._x = -_x;
    a._y = -_y;
    a._z = -_z;
  }

  /// @brief Normalize
  inline Vector3 &normalize() {
    double const length_ = get_length();
    assert(length_ != double(0));
    double const inv_length(double(1) / length_);
    _x *= inv_length;
    _y *= inv_length;
    _z *= inv_length;
    return *this;
  }

public: // Properties: accessors
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

public: // Properties: double assignment
  /// @brief x assignment
  inline void set_x(const double &x_a) { _x = x_a; }

  /// @brief y assignment
  inline void set_y(const double &y_a) { _y = y_a; }

  /// @brief z assignment
  inline void set_z(const double &z_a) { _z = z_a; }

private: // Methods
  /// @brief square( t ) == t * t
  inline static double square(const double &t) { return t * t; }

private: // Fields
  /// @brief Coordinates of the 3 coordinate vector
  double _x;
  double _y;
  double _z;

}; // Vector3

typedef std::vector<Vector3> Vector3s;

/// @brief generates a Vector3 object from a string with 3 numbers seperated by
/// by spaces.
inline Vector3 vector_from_str(const String &s) {
  Strings doubles = base::string::split(s, " ");
  Reals point;
  if (doubles.size() != 3) {
    String msg = "vector must contain 3 numbers. From string got " +
                 std::to_string(doubles.size());
    base::log_and_throw<base::InputException>(msg);
  }
  for (auto &i : doubles) {
    point.push_back(std::stod(i));
  }
  return Vector3{point};
}

inline void vectors_from_str(const String &s, Vector3s &vecs /* return */) {
  Strings doubles = base::string::split(s, " ");
  Reals point = {0, 0, 0};
  int pos = 0;
  for (auto const &d : doubles) {
    point[pos] = std::stod(d);
    pos += 1;
    if (pos == 3) {
      vecs.push_back(Vector3{point});
      pos = 0;
    }
  }
  // there are leftovers!
  if (pos != 0) {
    String msg = "vector must contain 3 numbers. From string got " +
                 std::to_string(doubles.size());
    base::log_and_throw<base::InputException>(msg);
  }
}

} // namespace math

#endif /* defined(__REDESIGNC__Vector3__) */