//
//  Matrix3x3.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 9/28/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__Matrix3x3__h
#define __RNAMake__Matrix3x3__h

#include <cmath>
#include <sstream>
#include <vector>
#include <iostream>

// RNAMake Headers
#include <base/string.hpp>
#include <base/types.hpp>
#include <math/transform.hpp>
#include <math/vector_3.hpp>
#include <base/exception.hpp>

namespace math {
// add comment here

class Matrix3x3 {


  // initiation ////////////////////////////////////////////////////////
public:
  inline Matrix3x3() = default;
  inline Matrix3x3(
          const double &xx, const double &xy, const double &xz,
          const double &yx, const double &yy, const double &yz,
          const double &zx, const double &zy, const double &zz)
      : _xx(xx), _xy(xy), _xz(xz),
        _yx(yx), _yy(yy), _yz(yz),
        _zx(zx), _zy(zy), _zz(zz) {}
  /// @brief Uniform value constructor
  inline explicit Matrix3x3(const double &t)
      : _xx(t), _xy(t), _xz(t),
        _yx(t), _yy(t), _yz(t),
        _zx(t), _zy(t), _zz(t) {
  }

  /// @brief Destructor
  inline ~Matrix3x3() = default;

public:
  inline String const get_str() const {
    return std::to_string(_xx) + " " + std::to_string(_xy) + " " + std::to_string(_xz) + " "
           + std::to_string(_yx) + " " + std::to_string(_yy) + " " + std::to_string(_yz) + " "
           + std::to_string(_zx) + " " + std::to_string(_zy) + " " + std::to_string(_zz);
  }

  [[nodiscard]] inline String const get_str_readable() const {
    auto s = String();
    s += "[" + std::to_string(_xx) + ", " + std::to_string(_xy) + ", " + std::to_string(_xz) + ",\n";
    s += " " + std::to_string(_yx) + ", " + std::to_string(_yy) + ", " + std::to_string(_yz) + ",\n";
    s += " " + std::to_string(_zx) + ", " + std::to_string(_zy) + ", " + std::to_string(_zz) + "]";
    return s;
  }

public:
  inline static Matrix3x3 identity() {
    return {1, 0, 0,
            0, 1, 0,
            0, 0, 1};
  }

  /// @brief Copy assignment
  inline Matrix3x3 &operator=(Matrix3x3 const &m) {
    _xx = m._xx;
    _xy = m._xy;
    _xz = m._xz;
    _yx = m._yx;
    _yy = m._yy;
    _yz = m._yz;
    _zx = m._zx;
    _zy = m._zy;
    _zz = m._zz;
    return *this;
  }

  /// @brief += Matrix3x3
  inline Matrix3x3 &operator+=(Matrix3x3 const &m) {
    _xx += m._xx; _xy += m._xy; _xz += m._xz;
    _yx += m._yx; _yy += m._yy; _yz += m._yz;
    _zx += m._zx; _zy += m._zy; _zz += m._zz;
    return *this;
  }

  inline Matrix3x3 &operator-=(Matrix3x3 const &m) {
    _xx -= m._xx; _xy -= m._xy; _xz -= m._xz;
    _yx -= m._yx; _yy -= m._yy; _yz -= m._yz;
    _zx -= m._zx; _zy -= m._zy; _zz -= m._zz;
    return *this;
  }

public:// Assignment: scalar
  /// @brief = Value
  inline Matrix3x3 &operator=(double const &t) {
    _xx = _xy = _xz = t;
    _yx = _yy = _yz = t;
    _zx = _zy = _zz = t;
    return *this;
  }

  /// @brief += Value
  inline Matrix3x3 &operator+=(Real const &t) {
    _xx += t; _xy += t; _xz += t;
    _yx += t; _yy += t; _yz += t;
    _zx += t; _zy += t; _zz += t;
    return *this;
  }

  /// @brief -= Value
  inline Matrix3x3 &operator-=(Real const &t) {
    _xx -= t; _xy -= t; _xz -= t;
    _yx -= t; _yy -= t; _yz -= t;
    _zx -= t; _zy -= t; _zz -= t;
    return *this;
  }

public: // Methods: basic mathematical
  /// @brief Matrix3x3 + Matrix3x3
  friend inline Matrix3x3 operator+(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(
            a._xx + b._xx, a._xy + b._xy, a._xz + b._xz,
            a._yx + b._yx, a._yy + b._yy, a._yz + b._yz,
            a._zx + b._zx, a._zy + b._zy, a._zz + b._zz);
  }

  /// @brief Matrix3x3 + Real
  friend inline Matrix3x3 operator+(Matrix3x3 const &m, Real const &t) {
    return Matrix3x3(
            m._xx + t, m._xy + t, m._xz + t,
            m._yx + t, m._yy + t, m._yz + t,
            m._zx + t, m._zy + t, m._zz + t);
  }

  /// @brief Real + Matrix3x3
  friend inline Matrix3x3 operator+(Real const &t, Matrix3x3 const &m) {
    return Matrix3x3(
            t + m._xx, t + m._xy, t + m._xz,
            t + m._yx, t + m._yy, t + m._yz,
            t + m._zx, t + m._zy, t + m._zz);
  }

  /// @brief Matrix3x3 - Matrix3x3
  friend inline Matrix3x3 operator-(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(
            a._xx - b._xx, a._xy - b._xy, a._xz - b._xz,
            a._yx - b._yx, a._yy - b._yy, a._yz - b._yz,
            a._zx - b._zx, a._zy - b._zy, a._zz - b._zz);
  }

  /// @brief Matrix3x3 - Value
  friend inline Matrix3x3 operator-(Matrix3x3 const &m, Real const &t) {
    return Matrix3x3(
            m._xx - t, m._xy - t, m._xz - t,
            m._yx - t, m._yy - t, m._yz - t,
            m._zx - t, m._zy - t, m._zz - t);
  }

  /// @brief Value - Matrix3x3
  friend inline Matrix3x3 operator-(Real const &t, Matrix3x3 const &m) {
    return Matrix3x3(
            t - m._xx, t - m._xy, t - m._xz,
            t - m._yx, t - m._yy, t - m._yz,
            t - m._zx, t - m._zy, t - m._zz);
  }

  /// @brief Matrix3x3 * Matrix3x3
  friend inline Matrix3x3 operator*(Matrix3x3 const &a, const Matrix3x3 &b) {
    return Matrix3x3(
            // First row
            (a._xx * b._xx) + (a._xy * b._yx) + (a._xz * b._zx),
            (a._xx * b._xy) + (a._xy * b._yy) + (a._xz * b._zy),
            (a._xx * b._xz) + (a._xy * b._yz) + (a._xz * b._zz),

            // Second row
            (a._yx * b._xx) + (a._yy * b._yx) + (a._yz * b._zx),
            (a._yx * b._xy) + (a._yy * b._yy) + (a._yz * b._zy),
            (a._yx * b._xz) + (a._yy * b._yz) + (a._yz * b._zz),

            // Third row
            (a._zx * b._xx) + (a._zy * b._yx) + (a._zz * b._zx),
            (a._zx * b._xy) + (a._zy * b._yy) + (a._zz * b._zy),
            (a._zx * b._xz) + (a._zy * b._yz) + (a._zz * b._zz));
  }

  /// @brief Matrix3x3 * Real
  friend inline Matrix3x3 operator*(Matrix3x3 const &m, Real const &t) {
    return Matrix3x3(
            m._xx * t, m._xy * t, m._xz * t,
            m._yx * t, m._yy * t, m._yz * t,
            m._zx * t, m._zy * t, m._zz * t);
  }

  /// @brief Value * Matrix3x3
  friend inline Matrix3x3 operator*(Real const &t, Matrix3x3 const &m) {
    return Matrix3x3(
            t * m._xx, t * m._xy, t * m._xz,
            t * m._yx, t * m._yy, t * m._yz,
            t * m._zx, t * m._zy, t * m._zz);
  }

  /// @brief Matrix3x3 / Value
  friend inline Matrix3x3 operator/(Matrix3x3 const &m, Real const &t) {
    if (t == 0) {
      String msg = "trying to divide a vector by zero!";
      base::log_and_throw<base::MathException>(msg);
    }
    Real const inv_t(Real(1) / t);
    return Matrix3x3(
            m._xx * inv_t, m._xy * inv_t, m._xz * inv_t,
            m._yx * inv_t, m._yy * inv_t, m._yz * inv_t,
            m._zx * inv_t, m._zy * inv_t, m._zz * inv_t);
  }

public:
  /// @brief Transpose
  inline Matrix3x3 &transpose() {
    Real temp = _xy;
    _xy = _yx;
    _yx = temp;

    temp = _xz;
    _xz = _zx;
    _zx = temp;

    temp = _yz;
    _yz = _zy;
    _zy = temp;

    return *this;
  }

  inline void transpose(Matrix3x3 &a) {
    _xx = a._xx;
    _yy = a._yy;
    _zz = a._zz;

    _yx = a._xy;
    _xy = a._yx;
    _zx = a._xz;
    _xz = a._zx;
    _yz = a._zy;
    _zy = a._yz;
  }

  // multiplies a matrix and vector and spits out the value
  inline Vector3 dot(const Vector3 &v) {
    Vector3 new_v = {0, 0, 0};
    new_v.set_x(_xx * v.get_x() + _xy * v.get_y() + _xz * v.get_z());
    new_v.set_y(_yx * v.get_x() + _yy * v.get_y() + _yz * v.get_z());
    new_v.set_z(_zx * v.get_x() + _zy * v.get_y() + _zz * v.get_z());
    return new_v;
  }

  // multiplies a matrix with a vector and saves the value
  inline void dot(Vector3 const &v, Vector3 &vr /* return */ ) {
    vr.set_x(get_xx() * v.get_x() + get_xy() * v.get_y() + get_xz() * v.get_z());
    vr.set_y(get_yx() * v.get_x() + get_yy() * v.get_y() + get_yz() * v.get_z());
    vr.set_z(get_zx() * v.get_x() + get_zy() * v.get_y() + get_zz() * v.get_z());
  }

  // multiplies two matrices and saves the value
  inline void dot(const Matrix3x3 &b, Matrix3x3 &c /* return */) {
    c.set_xx(_xx * b.get_xx() + _xy * b.get_yx() + _xz * b.get_zx());
    c.set_xy(_xx * b.get_xy() + _xy * b.get_yy() + _xz * b.get_zy());
    c.set_xz(_xx * b.get_xz() + _xy * b.get_yz() + _xz * b.get_zz());

    c.set_yx(_yx * b.get_xx() + _yy * b.get_yx() + _yz * b.get_zx());
    c.set_yy(_yx * b.get_xy() + _yy * b.get_yy() + _yz * b.get_zy());
    c.set_yz(_yx * b.get_xz() + _yy * b.get_yz() + _yz * b.get_zz());

    c.set_zx(_zx * b.get_xx() + _zy * b.get_yx() + _zz * b.get_zx());
    c.set_zy(_zx * b.get_xy() + _zy * b.get_yy() + _zz * b.get_zy());
    c.set_zz(_zx * b.get_xz() + _zy * b.get_yz() + _zz * b.get_zz());
  }

  // multiplies two matrices and spits out the value
  inline Matrix3x3 dot(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(
            a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() + a.get_xz() * b.get_zx(),
            a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() + a.get_xz() * b.get_zy(),
            a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() + a.get_xz() * b.get_zz(),

            a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() + a.get_yz() * b.get_zx(),
            a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() + a.get_yz() * b.get_zy(),
            a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() + a.get_yz() * b.get_zz(),

            a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() + a.get_zz() * b.get_zx(),
            a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() + a.get_zz() * b.get_zy(),
            a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() + a.get_zz() * b.get_zz());
  }

  inline double difference(Matrix3x3 const &b) const {
    double dist = 0.0f;
    dist += std::abs(_xx - b._xx);
    dist += std::abs(_xy - b._xy);
    dist += std::abs(_xz - b._xz);
    dist += std::abs(_yx - b._yx);
    dist += std::abs(_yy - b._yy);
    dist += std::abs(_yz - b._yz);
    dist += std::abs(_zx - b._zx);
    dist += std::abs(_zy - b._zy);
    dist += std::abs(_zz - b._zz);

    return dist;
  }

  inline Matrix3x3 get_flip_orientation() const {
    return Matrix3x3(
            _xx, _xy, _xz,
            -_yx, -_yy, -_yz,
            -_zx, -_zy, -_zz);
  }


  // TODO find out what unitarize is and write unittests
  inline Matrix3x3 get_unitarize() const {
    auto m = Matrix3x3(
            _xx, _xy, _xz,
            _yx, _yy, _yz,
            _zx, _zy, _zz
            );

    // R[0] /= math.sqrt(R[0].dot(R[0]))
    double dot = sqrt(_xx * _xx + _xy * _xy + _xz * _xz);
    m._xx /= dot;
    m._xy /= dot;
    m._xz /= dot;
    // R[1] -= R[1].dot(R[0]) * R[0]
    dot = _yx * m._xx + _yy * m._xy + _yz * m._xz;
    m._yx -= dot * m._xx;
    m._yy -= dot * m._xy;
    m._yz -= dot * m._xz;
    // R[1] /= math.sqrt(R[1].dot(R[1]))
    dot = sqrt(m._yx * m._yx + m._yy * m._yy + m._yz * m._yz);
    m._yx /= dot;
    m._yy /= dot;
    m._yz /= dot;
    // R[2] -= R[2].dot(R[0]) * R[0]
    dot = m._zx * m._xx + m._zy * m._xy + m._zz * m._xz;
    m._zx -= dot * m._xx;
    m._zy -= dot * m._xy;
    m._zz -= dot * m._xz;
    // R[2] -= R[2].dot(R[1]) * R[1]
    dot = m._zx * m._yx + m._zy * m._yy + m._zz * m._yz;
    m._zx -= dot * m._yx;
    m._zy -= dot * m._yy;
    m._zz -= dot * m._yz;
    // R[2] /= math.sqrt(R[2].dot(R[2]))
    dot = sqrt(m._zx * m._zx + m._zy * m._zy + m._zz * m._zz);
    m._zx /= dot;
    m._zy /= dot;
    m._zz /= dot;

    return m;
  }
// TODO find out what it does and write unittests
  inline void unitarize() {
    double dot = sqrt(_xx * _xx + _xy * _xy + _xz * _xz);
    _xx /= dot;
    _xy /= dot;
    _xz /= dot;
    // R[1] -= R[1].dot(R[0]) * R[0]
    dot = _yx * _xx + _yy * _xy + _yz * _xz;
    _yx -= dot * _xx;
    _yy -= dot * _xy;
    _yz -= dot * _xz;
    // R[1] /= math.sqrt(R[1].dot(R[1]))
    dot = sqrt(_yx * _yx + _yy * _yy + _yz * _yz);
    _yx /= dot;
    _yy /= dot;
    _yz /= dot;
    // R[2] -= R[2].dot(R[0]) * R[0]
    dot = _zx * _xx + _zy * _xy + _zz * _xz;
    _zx -= dot * _xx;
    _zy -= dot * _xy;
    _zz -= dot * _xz;
    // R[2] -= R[2].dot(R[1]) * R[1]
    dot = _zx * _yx + _zy * _yy + _zz * _yz;
    _zx -= dot * _yx;
    _zy -= dot * _yy;
    _zz -= dot * _yz;
    // R[2] /= math.sqrt(R[2].dot(R[2]))
    dot = sqrt(_zx * _zx + _zy * _zy + _zz * _zz);
    _zx /= dot;
    _zy /= dot;
    _zz /= dot;
  }

public:// Properties: scalars
  /// @brief Value xx const
  [[nodiscard]] inline Real const &get_xx() const { return _xx; }

  /// @brief Value xy const
  [[nodiscard]] inline Real const &get_xy() const { return _xy; }

  /// @brief Value xz const
  [[nodiscard]] inline Real const &get_xz() const { return _xz; }

  /// @brief Value yx const
  [[nodiscard]] inline Real const &get_yx() const { return _yx; }

  /// @brief Value yy const
  [[nodiscard]] inline Real const &get_yy() const { return _yy; }

  /// @brief Value yz const
  [[nodiscard]] inline Real const &get_yz() const { return _yz; }

  /// @brief Value zx const
  [[nodiscard]] inline Real const &get_zx() const { return _zx; }

  /// @brief Value zy const
  [[nodiscard]] inline Real const &get_zy() const { return _zy; }

  /// @brief Value zz const
  [[nodiscard]] inline Real const &get_zz() const { return _zz; }

public:// Properties: value assignment
  /// @brief xx assignment
  inline void set_xx(Real const &_xxa) { _xx = _xxa; }

  /// @brief xy assignment
  inline void set_xy(Real const &_xya) { _xy = _xya; }

  /// @brief xz assignment
  inline void set_xz(Real const &_xza) { _xz = _xza; }

  /// @brief yx assignment
  inline void set_yx(Real const &_yxa) { _yx = _yxa; }

  /// @brief yy assignment
  inline void set_yy(Real const &_yya) { _yy = _yya; }

  /// @brief yz assignment
  inline void set_yz(Real const &_yza) { _yz = _yza; }

  /// @brief zx assignment
  inline void set_zx(Real const &_zxa) { _zx = _zxa; }

  /// @brief zy assignment
  inline void set_zy(Real const &_zya) { _zy = _zya; }

  /// @brief zz assignment
  inline void set_zz(Real const &_zza) { _zz = _zza; }

  inline Matrix3x3 get_transposed() const {
    return Matrix3x3(
            _xx, _yx, _zx,
            _xy, _yy, _zy,
            _xz, _yz, _zz);
  }

  typedef std::vector<Matrix3x3> Matrix3x3s;

  // TODO Remove old code and check for duplicates

  //
  /*
  inline void dot_vectors(Vector3s const &v, Vector3s &vr) {
    int i;
    for (i = 0; i < v.size(); i++) {
      dot_vector(v[i], vr[i]);
    }
  }
   */

  inline Matrix3x3 transform_1() {
    return Matrix3x3(
            get_xx(), get_xy(), get_xz(),
            -get_yx(), -get_yy(), -get_yz(),
            -get_zx(), -get_zy(), -get_zz());
  }


  template<typename T>
  std::ostream &operator<<(std::ostream &stream) {
    stream << "(" << get_xx() << ", " << get_xy() << ", " << get_xz() << ")" << std::endl;
    stream << "(" << get_yx() << ", " << get_yy() << ", " << get_yz() << ")" << std::endl;
    stream << "(" << get_zx() << ", " << get_zy() << ", " << get_zz() << ")" << std::endl;
    return stream;
  }

  inline String matrix_to_str(Matrix3x3 const &m) {
    std::stringstream ss;
    ss << m.get_xx() << " " << m.get_xy() << " " << m.get_xz() << " ";
    ss << m.get_yx() << " " << m.get_yy() << " " << m.get_yz() << " ";
    ss << m.get_zx() << " " << m.get_zy() << " " << m.get_zz() << " ";
    return ss.str();
  }


private:
  Real _xx, _xy, _xz;
  Real _yx, _yy, _yz;
  Real _zx, _zy, _zz;
};

Matrix3x3 matrix_from_str(const String &s);

// add end comment on the line above
} // namespace math

#endif