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
#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>

namespace math {

class Transform {
  friend class Matrix3x3;

public:
   inline void dot(Matrix3x3 const &b, Transform &c) const {
    c._xx = _xx * b.get_xx() + _xy * b.get_yx() + _xz * b.get_zx();
    c._xy = _xx * b.get_xy() + _xy * b.get_yy() + _xz * b.get_zy();
    c._xz = _xx * b.get_xz() + _xy * b.get_yz() + _xz * b.get_zz();

    c._yx = _yx * b.get_xx() + _yy * b.get_yx() + _yz * b.get_zx();
    c._yy = _yx * b.get_xy() + _yy * b.get_yy() + _yz * b.get_zy();
    c._yz = _yx * b.get_xz() + _yy * b.get_yz() + _yz * b.get_zz();

    c._zx = _zx * b.get_xx() + _zy * b.get_yx() + _zz * b.get_zx();
    c._zy = _zx * b.get_xy() + _zy * b.get_yy() + _zz * b.get_zy();
    c._zz = _zx * b.get_xz() + _zy * b.get_yz() + _zz * b.get_zz();
  }

 public:  // initialization ////////////////////////////////////////////////////
  /// @brief Default constructor
  inline Transform()
      : _xx(1.0), _yx(0.0), _zx(0.0), _px(0.0),
        _xy(0.0), _yy(1.0), _zy(0.0), _py(0.0),
        _xz(0.0), _yz(0.0), _zz(1.0), _pz(0.0) {}

  inline Transform(
          const double &xx, const double &xy, const double &xz,
          const double &yx, const double &yy, const double &yz,
          const double &zx, const double &zy, const double &zz,
          const double &px, const double &py, const double &pz)
      : _xx(xx), _xy(xy), _xz(xz),
        _yx(yx), _yy(yy), _yz(yz),
        _zx(zx), _zy(zy), _zz(zz),
        _px(px), _py(py), _pz(pz) {}

  inline Transform(Matrix3x3 const& r, Vector3 const& t) {
    rotation(r);
    translation(t);
  }

  /// @brief - destructor
  inline ~Transform() = default;

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

  [[nodiscard]] Vector3 xaxis() const { return Vector3(_xx, _xy, _xz); }

  [[nodiscard]] Vector3 yaxis() const { return Vector3(_yx, _yy, _yz); }

  [[nodiscard]] Vector3 zaxis() const { return Vector3(_zx, _zy, _zz); }

  [[nodiscard]] Vector3 translation() const { return Vector3(_px, _py, _pz); }

 public:
  [[nodiscard]] inline Matrix3x3 rotation() const {
    return Matrix3x3(
            _xx, _xy, _xz,
            _yx, _yy, _yz,
            _zx, _zy, _zz);
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

}  // namespace math

#endif /* defined(__REDESIGNC__Transform__) */
