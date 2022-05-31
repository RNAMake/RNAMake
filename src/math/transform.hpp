//
//  Transform.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/1/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Transform__
#define __REDESIGNC__Transform__

#include <iostream>

// RNAMake Headers

#include "matrix_3x3.hpp"
#include "vector_3.hpp"
#include <../math/matrix_3x3.hpp>
#include <../math/vector_3.hpp>

namespace math {
/*
class Transform {

 public:
  friend inline void dot(math::Matrix3x3 const& a, math::Matrix3x3 const& b, Transform& c) {
    c._xx = a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() + a.get_xz() * b.get_zx();
    c._xy = a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() + a.get_xz() * b.get_zy();
    c._xz = a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() + a.get_xz() * b.get_zz();

    c._yx = a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() + a.get_yz() * b.get_zx();
    c._yy = a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() + a.get_yz() * b.get_zy();
    c._yz = a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() + a.get_yz() * b.get_zz();

    c._zx = a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() + a.get_zz() * b.get_zx();
    c._zy = a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() + a.get_zz() * b.get_zy();
    c._zz = a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() + a.get_zz() * b.get_zz();
  }

 public:  // creation
  /// @brief Default constructor
  inline Transform()
      : _xx(1.0),
        _yx(0.0),
        _zx(0.0),
        _px(0.0),
        _xy(0.0),
        _yy(1.0),
        _zy(0.0),
        _py(0.0),
        _xz(0.0),
        _yz(0.0),
        _zz(1.0),
        _pz(0.0) {}

  inline Transform(math::Matrix3x3 const& r, Vector3 const& t) {
    rotation(r);
    translation(t);
  }

 public:
  // Accessors

  [[nodiscard]] float get_xx() const { return _xx; }

  [[nodiscard]] float get_xy() const { return _xy; }

  [[nodiscard]] float get_xz() const { return _xz; }

  [[nodiscard]] float get_yx() const { return _yx; }

  [[nodiscard]] float get_yy() const { return _yy; }

  [[nodiscard]] float get_yz() const { return _yz; }

  [[nodiscard]] float get_zx() const { return _zx; }

  [[nodiscard]] float get_zy() const { return _zy; }

  [[nodiscard]] float get_zz() const { return _zz; }

  [[nodiscard]] float get_px() const { return _px; }

  [[nodiscard]] float get_py() const { return _py; }

  [[nodiscard]] float get_pz() const { return _pz; }

  float get_xx() { return _xx; }

  float get_xy() { return _xy; }

  float get_xz() { return _xz; }

  float get_yx() { return _yx; }

  float get_yy() { return _yy; }

  float get_yz() { return _yz; }

  float get_zx() { return _zx; }

  float get_zy() { return _zy; }

  float get_zz() { return _zz; }

  float get_px() { return _px; }

  float get_py() { return _py; }

  float get_pz() { return _pz; }

  Vector3 xaxis() const { return Vector3(_xx, _xy, _xz); }

  Vector3 yaxis() const { return Vector3(_yx, _yy, _yz); }

  Vector3 zaxis() const { return Vector3(_zx, _zy, _zz); }

  Vector3 translation() const { return Vector3(_px, _py, _pz); }

 public:
  inline Matrix3x3 rotation() const {
    return Matrix3x3(_xx, _xy, _xz, _yx, _yy, _yz, _zx, _zy, _zz);
  }

  inline void rotation(Matrix3x3 const& m) {
    _xx = m.get_xx();
    _xy = m.get_xy();
    _xz = m.get_xz();
    _yx = m.get_yx();
    _yy = m.get_yy();
    _yz = m.get_yz();
    _zx = m.get_zx();
    _zy = m.get_zy();
    _zz = m.get_zz();
  }

  inline void translation(Vector3 const& v) {
    _px = v.get_x();
    _py = v.get_y();
    _pz = v.get_z();
  }

 private:
  float _xx, _yx, _zx, _px;
  float _xy, _yy, _zy, _py;
  float _xz, _yz, _zz, _pz;
};
*/
}  // namespace math

#endif /* defined(__REDESIGNC__Transform__) */
