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
#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>

namespace math {
/*
class Transform {
 public:

 public:
  friend inline void dot(math::Matrix3x3 const& a, math::Matrix3x3 const& b, Transform& c) {
    c._xx = a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() +
            a.get_xz() * b.get_zx();
    c._xy = a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() +
            a.get_xz() * b.get_zy();
    c._xz = a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() +
            a.get_xz() * b.get_zz();
    c._yx = a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() +
            a.get_yz() * b.get_zx();
    c._yy = a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() +
            a.get_yz() * b.get_zy();
    c._yz = a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() +
            a.get_yz() * b.get_zz();
    c._zx = a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() +
            a.get_zz() * b.get_zx();
    c._zy = a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() +
            a.get_zz() * b.get_zy();
    c._zz = a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() +
            a.get_zz() * b.get_zz();
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

  inline Transform(Matrix3x3 const& r, Point const& t) {
    rotation(r);
    translation(t);
  }

 public:
  // Accessors

  float get_xx() const { return _xx; }

  float xy() const { return _xy; }

  float xz() const { return _xz; }

  float yx() const { return _yx; }

  float yy() const { return _yy; }

  float yz() const { return _yz; }

  float zx() const { return _zx; }

  float zy() const { return _zy; }

  float zz() const { return _zz; }

  float px() const { return _px; }

  float py() const { return _py; }

  float pz() const { return _pz; }

  float get_xx() { return _xx; }

  float xy() { return _xy; }

  float xz() { return _xz; }

  float yx() { return _yx; }

  float yy() { return _yy; }

  float yz() { return _yz; }

  float zx() { return _zx; }

  float zy() { return _zy; }

  float zz() { return _zz; }

  float px() { return _px; }

  float py() { return _py; }

  float pz() { return _pz; }

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
