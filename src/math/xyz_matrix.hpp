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

// RNAMake Headers
#include <base/string.hpp>
#include <base/types.hpp>
#include <math/vector_3.hpp>

namespace math {
/*
class Matrix3x3 {
public:
  friend class Transform;

  friend inline void transpose(Matrix3x3 const &a, Matrix3x3 &b) {
    b._xx = a._xx;
    b._yy = a._yy;
    b.zz_ = a.zz_;

    b._yx = a._xy;
    b._xy = a._yx;
    b.zx_ = a.xz_;
    b.xz_ = a.zx_;
    b.yz_ = a.zy_;
    b.zy_ = a.yz_;
  }

  friend inline Vector3 dot_vector(Matrix3x3 const &m, Vector3 &v) {
    auto new_v = Vector3(0, 0, 0);
    new_v.set_x(m._xx * v.get_x() + m._yx * v.get_y() + m.zx_ * v.get_z());
    new_v.set_y(m._xy * v.get_x() + m._yy * v.get_y() + m.zy_ * v.get_z());
    new_v.set_z(m.xz_ * v.get_x() + m.yz_ * v.get_y() + m.zz_ * v.get_z());
    return new_v;
  }

public: // initiation ////////////////////////////////////////////////////////
  inline Matrix3x3() = default;
  inline Matrix3x3(Matrix3x3 const &m) = default;
  inline Matrix3x3(const double &xx, const double xy, const double xz,
                   const double yx, const double yy, const double yz,
                   const double zx, const double zy, const double zz)
      : _xx(xx), _xy(xy), xz_(xz), _yx(yx), _yy(yy), yz_(yz), zx_(zx), zy_(zy),
        zz_(zz) {}
  /// @brief Uniform value constructor
  inline explicit Matrix3x3(const double &t)
      : _xx(t), _xy(t), xz_(t), _yx(t), _yy(t), yz_(t), zx_(t), zy_(t), zz_(t) {
  }



  /// @brief Destructor
  inline ~Matrix3x3() = default;

public:
  inline String const get_str() const {
    return std::to_string(_xx) + " " + std::to_string(_xy) + " " +
           std::to_string(xz_) + " " + std::to_string(_yx) + " " +
           std::to_string(_yy) + " " + std::to_string(yz_) + " " +
           std::to_string(zx_) + " " + std::to_string(zy_) + " " +
           std::to_string(zz_);
  }

  [[nodiscard]] inline String const get_str_readable() const {
    auto s = String();
    s += "[" + std::to_string(_xx) + ", " + std::to_string(_xy) + ", " +
         std::to_string(xz_) + ",\n";
    s += " " + std::to_string(_yx) + ", " + std::to_string(_yy) + ", " +
         std::to_string(yz_) + ",\n";
    s += " " + std::to_string(zx_) + ", " + std::to_string(zy_) + ", " +
         std::to_string(zz_) + "]";
    return s;
  }

public:
  inline static Matrix3x3 identity() {
    return {1, 0, 0, 0, 1, 0, 0, 0, 1};
  }

  /// @brief Copy assignment
  inline Matrix3x3 &operator=(Matrix3x3 const &m) {
    _xx = m._xx;
    _xy = m._xy;
    xz_ = m.xz_;
    _yx = m._yx;
    _yy = m._yy;
    yz_ = m.yz_;
    zx_ = m.zx_;
    zy_ = m.zy_;
    zz_ = m.zz_;
    return *this;
  }

  /// @brief += Matrix3x3
  inline Matrix3x3 &operator+=(Matrix3x3 const &m) {
    _xx += m._xx;
    _xy += m._xy;
    xz_ += m.xz_;
    _yx += m._yx;
    _yy += m._yy;
    yz_ += m.yz_;
    zx_ += m.zx_;
    zy_ += m.zy_;
    zz_ += m.zz_;
    return *this;
  }

  inline Matrix3x3 &operator-=(Matrix3x3 const &m) {
    _xx -= m._xx;
    _xy -= m._xy;
    xz_ -= m.xz_;
    _yx -= m._yx;
    _yy -= m._yy;
    yz_ -= m.yz_;
    zx_ -= m.zx_;
    zy_ -= m.zy_;
    zz_ -= m.zz_;
    return *this;
  }

public: // Assignment: scalar
  /// @brief = Value
  inline Matrix3x3 &operator=(double const &t) {
    _xx = _xy = xz_ = t;
    _yx = _yy = yz_ = t;
    zx_ = zy_ = zz_ = t;
    return *this;
  }

  /// @brief += Value
  inline Matrix3x3 &operator+=(Value const &t) {
    _xx += t;
    _xy += t;
    xz_ += t;
    _yx += t;
    _yy += t;
    yz_ += t;
    zx_ += t;
    zy_ += t;
    zz_ += t;
    return *this;
  }

  /// @brief -= Value
  inline Matrix3x3 &operator-=(Value const &t) {
    _xx -= t;
    _xy -= t;
    xz_ -= t;
    _yx -= t;
    _yy -= t;
    yz_ -= t;
    zx_ -= t;
    zy_ -= t;
    zz_ -= t;
    return *this;
  }

public: // Methods: basic mathematical
  /// @brief Matrix3x3 + Matrix3x3
  friend inline Matrix3x3 operator+(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(a._xx + b._xx, a._xy + b._xy, a.xz_ + b.xz_, a._yx + b._yx,
                     a._yy + b._yy, a.yz_ + b.yz_, a.zx_ + b.zx_, a.zy_ + b.zy_,
                     a.zz_ + b.zz_);
  }

  /// @brief Matrix3x3 + Value
  friend inline Matrix3x3 operator+(Matrix3x3 const &m, Value const &t) {
    return Matrix3x3(m._xx + t, m._xy + t, m.xz_ + t, m._yx + t, m._yy + t,
                     m.yz_ + t, m.zx_ + t, m.zy_ + t, m.zz_ + t);
  }

  /// @brief Value + Matrix3x3
  friend inline Matrix3x3 operator+(Value const &t, Matrix3x3 const &m) {
    return Matrix3x3(t + m._xx, t + m._xy, t + m.xz_, t + m._yx, t + m._yy,
                     t + m.yz_, t + m.zx_, t + m.zy_, t + m.zz_);
  }

  /// @brief Matrix3x3 - Matrix3x3
  friend inline Matrix3x3 operator-(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(a._xx - b._xx, a._xy - b._xy, a.xz_ - b.xz_, a._yx - b._yx,
                     a._yy - b._yy, a.yz_ - b.yz_, a.zx_ - b.zx_, a.zy_ - b.zy_,
                     a.zz_ - b.zz_);
  }

  /// @brief Matrix3x3 - Value
  friend inline Matrix3x3 operator-(Matrix3x3 const &m, Value const &t) {
    return Matrix3x3(m._xx - t, m._xy - t, m.xz_ - t, m._yx - t, m._yy - t,
                     m.yz_ - t, m.zx_ - t, m.zy_ - t, m.zz_ - t);
  }

  /// @brief Value - Matrix3x3
  friend inline Matrix3x3 operator-(Value const &t, Matrix3x3 const &m) {
    return Matrix3x3(t - m._xx, t - m._xy, t - m.xz_, t - m._yx, t - m._yy,
                     t - m.yz_, t - m.zx_, t - m.zy_, t - m.zz_);
  }

  /// @brief Matrix3x3 * Matrix3x3
  friend inline Matrix3x3 operator*(Matrix3x3 const &a, Matrix3x3 const &b) {
    return Matrix3x3(
        // First row
        (a._xx * b._xx) + (a._xy * b._yx) + (a.xz_ * b.zx_),
        (a._xx * b._xy) + (a._xy * b._yy) + (a.xz_ * b.zy_),
        (a._xx * b.xz_) + (a._xy * b.yz_) + (a.xz_ * b.zz_),

        // Second row
        (a._yx * b._xx) + (a._yy * b._yx) + (a.yz_ * b.zx_),
        (a._yx * b._xy) + (a._yy * b._yy) + (a.yz_ * b.zy_),
        (a._yx * b.xz_) + (a._yy * b.yz_) + (a.yz_ * b.zz_),

        // Third row
        (a.zx_ * b._xx) + (a.zy_ * b._yx) + (a.zz_ * b.zx_),
        (a.zx_ * b._xy) + (a.zy_ * b._yy) + (a.zz_ * b.zy_),
        (a.zx_ * b.xz_) + (a.zy_ * b.yz_) + (a.zz_ * b.zz_));
  }

  /// @brief Matrix3x3 * Value
  friend inline Matrix3x3 operator*(Matrix3x3 const &m, Value const &t) {
    return Matrix3x3(m._xx * t, m._xy * t, m.xz_ * t, m._yx * t, m._yy * t,
                     m.yz_ * t, m.zx_ * t, m.zy_ * t, m.zz_ * t);
  }

  /// @brief Value * Matrix3x3
  friend inline Matrix3x3 operator*(Value const &t, Matrix3x3 const &m) {
    return Matrix3x3(t * m._xx, t * m._xy, t * m.xz_, t * m._yx, t * m._yy,
                     t * m.yz_, t * m.zx_, t * m.zy_, t * m.zz_);
  }

  /// @brief Matrix3x3 / Value
  friend inline Matrix3x3 operator/(Matrix3x3 const &m, Value const &t) {
    assert(t != Value(0));
    Value const inv_t(Value(1) / t);
    return Matrix3x3(m._xx * inv_t, m._xy * inv_t, m.xz_ * inv_t, m._yx * inv_t,
                     m._yy * inv_t, m.yz_ * inv_t, m.zx_ * inv_t, m.zy_ * inv_t,
                     m.zz_ * inv_t);
  }

public:
  /// @brief Transpose
  inline Matrix3x3 &transpose() {
    Value temp = _xy;
    _xy = _yx;
    _yx = temp;

    temp = xz_;
    xz_ = zx_;
    zx_ = temp;

    temp = yz_;
    yz_ = zy_;
    zy_ = temp;

    return *this;
  }

  inline double difference(Matrix3x3<T> const &b) const {
    double dist = 0.0f;
    dist += std::abs(_xx - b._xx);
    dist += std::abs(_xy - b._xy);
    dist += std::abs(xz_ - b.xz_);
    dist += std::abs(_yx - b._yx);
    dist += std::abs(_yy - b._yy);
    dist += std::abs(yz_ - b.yz_);
    dist += std::abs(zx_ - b.zx_);
    dist += std::abs(zy_ - b.zy_);
    dist += std::abs(zz_ - b.zz_);

    return dist;
  }

  inline Matrix3x3 get_flip_orientation() const {
    return Matrix3x3(_xx, _xy, xz_, -_yx, -_yy, -yz_, -zx_, -zy_, -zz_);
  }

  inline Matrix3x3 get_unitarize() const {
    auto m = Matrix3x3(_xx, _xy, xz_, _yx, _yy, yz_, zx_, zy_, zz_);

    // R[0] /= math.sqrt(R[0].dot(R[0]))
    double dot = sqrt(_xx * _xx + _xy * _xy + xz_ * xz_);
    m._xx /= dot;
    m._xy /= dot;
    m.xz_ /= dot;
    // R[1] -= R[1].dot(R[0]) * R[0]
    dot = _yx * m._xx + _yy * m._xy + yz_ * m.xz_;
    m._yx -= dot * m._xx;
    m._yy -= dot * m._xy;
    m.yz_ -= dot * m.xz_;
    // R[1] /= math.sqrt(R[1].dot(R[1]))
    dot = sqrt(m._yx * m._yx + m._yy * m._yy + m.yz_ * m.yz_);
    m._yx /= dot;
    m._yy /= dot;
    m.yz_ /= dot;
    // R[2] -= R[2].dot(R[0]) * R[0]
    dot = m.zx_ * m._xx + m.zy_ * m._xy + m.zz_ * m.xz_;
    m.zx_ -= dot * m._xx;
    m.zy_ -= dot * m._xy;
    m.zz_ -= dot * m.xz_;
    // R[2] -= R[2].dot(R[1]) * R[1]
    dot = m.zx_ * m._yx + m.zy_ * m._yy + m.zz_ * m.yz_;
    m.zx_ -= dot * m._yx;
    m.zy_ -= dot * m._yy;
    m.zz_ -= dot * m.yz_;
    // R[2] /= math.sqrt(R[2].dot(R[2]))
    dot = sqrt(m.zx_ * m.zx_ + m.zy_ * m.zy_ + m.zz_ * m.zz_);
    m.zx_ /= dot;
    m.zy_ /= dot;
    m.zz_ /= dot;

    return m;
  }

  inline void unitarize() {
    double dot = sqrt(_xx * _xx + _xy * _xy + xz_ * xz_);
    _xx /= dot;
    _xy /= dot;
    xz_ /= dot;
    // R[1] -= R[1].dot(R[0]) * R[0]
    dot = _yx * _xx + _yy * _xy + yz_ * xz_;
    _yx -= dot * _xx;
    _yy -= dot * _xy;
    yz_ -= dot * xz_;
    // R[1] /= math.sqrt(R[1].dot(R[1]))
    dot = sqrt(_yx * _yx + _yy * _yy + yz_ * yz_);
    _yx /= dot;
    _yy /= dot;
    yz_ /= dot;
    // R[2] -= R[2].dot(R[0]) * R[0]
    dot = zx_ * _xx + zy_ * _xy + zz_ * xz_;
    zx_ -= dot * _xx;
    zy_ -= dot * _xy;
    zz_ -= dot * xz_;
    // R[2] -= R[2].dot(R[1]) * R[1]
    dot = zx_ * _yx + zy_ * _yy + zz_ * yz_;
    zx_ -= dot * _yx;
    zy_ -= dot * _yy;
    zz_ -= dot * yz_;
    // R[2] /= math.sqrt(R[2].dot(R[2]))
    dot = sqrt(zx_ * zx_ + zy_ * zy_ + zz_ * zz_);
    zx_ /= dot;
    zy_ /= dot;
    zz_ /= dot;
  }

public: // Properties: scalars
  /// @brief Value xx const
  inline Value const &get_xx() const { return _xx; }

  /// @brief Value xy const
  inline Value const &get_xy() const { return _xy; }

  /// @brief Value xz const
  inline Value const &get_xz() const { return xz_; }

  /// @brief Value yx const
  inline Value const &get_yx() const { return _yx; }

  /// @brief Value yy const
  inline Value const &get_yy() const { return _yy; }

  /// @brief Value yz const
  inline Value const &get_yz() const { return yz_; }

  /// @brief Value zx const
  inline Value const &get_zx() const { return zx_; }

  /// @brief Value zy const
  inline Value const &get_zy() const { return zy_; }

  /// @brief Value zz const
  inline Value const &get_zz() const { return zz_; }

public: // Properties: value assignment
  /// @brief xx assignment
  inline void set_xx(Value const &_xxa) { _xx = _xxa; }

  /// @brief xy assignment
  inline void set_xy(Value const &_xya) { _xy = _xya; }

  /// @brief xz assignment
  inline void set_xz(Value const &xz_a) { xz_ = xz_a; }

  /// @brief yx assignment
  inline void set_yx(Value const &_yxa) { _yx = _yxa; }

  /// @brief yy assignment
  inline void set_yy(Value const &_yya) { _yy = _yya; }

  /// @brief yz assignment
  inline void set_yz(Value const &yz_a) { yz_ = yz_a; }

  /// @brief zx assignment
  inline void set_zx(Value const &zx_a) { zx_ = zx_a; }

  /// @brief zy assignment
  inline void set_zy(Value const &zy_a) { zy_ = zy_a; }

  /// @brief zz assignment
  inline void set_zz(Value const &zz_a) { zz_ = zz_a; }

  inline Matrix3x3<T> get_transposed() const {
    return Matrix3x3(_xx, _yx, zx_, _xy, _yy, zy_, xz_, yz_, zz_);
  }

private:
  Value _xx, _xy, xz_;
  Value _yx, _yy, yz_;
  Value zx_, zy_, zz_;
};

typedef Matrix3x3<double> Matrix;
typedef std::vector<Matrix> Matrices;

// TODO Remove old code and ask joe about transform_1
template <typename T>
inline Vector3 operator*(Matrix3x3<T> const &m, Vector3 const &v) {
  return Vector3(
      m.xx() * v.get_x() + m.get_xy() * v.get_y() + m.get_xz() * v.get_z(),
      m.get_yx() * v.get_x() + m.get_yy() * v.get_y() + m.get_yz() * v.get_z(),
      m.get_zx() * v.get_x() + m.get_zy() * v.get_y() + m.get_zz() * v.get_z());
}

inline void dot_vector(Matrix const &m, Vector3 const &v, Vector3 &vr) {
  vr.set_x(m.get_xx() * v.get_x() + m.get_yx() * v.get_y() +
           m.get_zx() * v.get_z());
  vr.set_y(m.get_xy() * v.get_x() + m.get_yy() * v.get_y() +
           m.get_zy() * v.get_z());
  vr.set_z(m.get_xz() * v.get_x() + m.get_yz() * v.get_y() +
           m.get_zz() * v.get_z());
}

inline Vector3 dot_vector(Matrix const &m, Vector3 const &v) {
  auto vr = Vector3(0, 0, 0);
  vr.set_x(m.get_xx() * v.get_x() + m.get_yx() * v.get_y() +
           m.get_zx() * v.get_z());
  vr.set_y(m.get_xy() * v.get_x() + m.get_yy() * v.get_y() +
           m.get_zy() * v.get_z());
  vr.set_z(m.get_xz() * v.get_x() + m.get_yz() * v.get_y() +
           m.get_zz() * v.get_z());
  return vr;
}

inline void dot_vectors(Matrix const &m, Vector3s const &v, Vector3s &vr) {
  int i;
  for (i = 0; i < v.size(); i++) {
    dot_vector(m, v[i], vr[i]);
  }
}

template <typename T> inline Matrix3x3<T> transform_1(Matrix3x3<T> const &m) {
  return Matrix3x3<T>(m.get_xx(), m.get_xy(), m.get_xz(), -m.get_yx(),
                      -m.get_yy(), -m.get_yz(), -m.get_zx(), -m.get_zy(),
                      -m.get_zz());
}

inline void dot(Matrix const &a, Matrix const &b, Matrix &c) {
  c.set_xx(a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() +
           a.get_xz() * b.get_zx());
  c.set_xy(a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() +
           a.get_xz() * b.get_zy());
  c.set_xz(a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() +
           a.get_xz() * b.get_zz());
  c.set_yx(a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() +
           a.get_yz() * b.get_zx());
  c.set_yy(a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() +
           a.get_yz() * b.get_zy());
  c.set_yz(a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() +
           a.get_yz() * b.get_zz());
  c.set_zx(a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() +
           a.get_zz() * b.get_zx());
  c.set_zy(a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() +
           a.get_zz() * b.get_zy());
  c.set_zz(a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() +
           a.get_zz() * b.get_zz());
}

inline Matrix dot(Matrix const &a, Matrix const &b) {
  return Matrix(a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() +
                    a.get_xz() * b.get_zx(),
                a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() +
                    a.get_xz() * b.get_zy(),
                a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() +
                    a.get_xz() * b.get_zz(),
                a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() +
                    a.get_yz() * b.get_zx(),
                a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() +
                    a.get_yz() * b.get_zy(),
                a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() +
                    a.get_yz() * b.get_zz(),
                a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() +
                    a.get_zz() * b.get_zx(),
                a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() +
                    a.get_zz() * b.get_zy(),
                a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() +
                    a.get_zz() * b.get_zz());
}

template <typename T>
std::ostream &operator<<(std::ostream &stream, Matrix3x3<T> const &v) {
  stream << "(" << v.get_xx() << ", " << v.get_xy() << ", " << v.get_xz() << ")"
         << std::endl;
  stream << "(" << v.get_yx() << ", " << v.get_yy() << ", " << v.get_yz() << ")"
         << std::endl;
  stream << "(" << v.get_zx() << ", " << v.get_zy() << ", " << v.get_zz() << ")"
         << std::endl;
  return stream;
}

inline String matrix_to_str(Matrix const &m) {
  std::stringstream ss;
  ss << m.get_xx() << " " << m.get_xy() << " " << m.get_xz() << " ";
  ss << m.get_yx() << " " << m.get_yy() << " " << m.get_yz() << " ";
  ss << m.get_zx() << " " << m.get_zy() << " " << m.get_zz() << " ";
  return ss.str();
}

  inline Matrix3x3(String const &s) {
  auto v = base::string::split(s, " ");
  assert(v.size() > 8);
  _xx = std::stod(v[0]);
  _xy = std::stod(v[1]);
  xz_ = std::stod(v[2]);
  _yx = std::stod(v[3]);
  _yy = std::stod(v[4]);
  yz_ = std::stod(v[5]);
  zx_ = std::stod(v[6]);
  zy_ = std::stod(v[7]);
  zz_ = std::stod(v[8]);
}

*/
} // namespace math

#endif