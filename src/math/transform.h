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
#include "math/xyz_matrix.h"
#include "math/xyz_vector.h"

namespace math {

class Transform {
 public:
  template <typename>
  friend class xyzMatrix;

  template <typename>
  friend class xyzVector;

 public:
  friend inline void dot(Matrix const& a, Matrix const& b, Transform& c) {
    c.xx_ = a.get_xx() * b.get_xx() + a.get_xy() * b.get_yx() +
            a.get_xz() * b.get_zx();
    c.xy_ = a.get_xx() * b.get_xy() + a.get_xy() * b.get_yy() +
            a.get_xz() * b.get_zy();
    c.xz_ = a.get_xx() * b.get_xz() + a.get_xy() * b.get_yz() +
            a.get_xz() * b.get_zz();
    c.yx_ = a.get_yx() * b.get_xx() + a.get_yy() * b.get_yx() +
            a.get_yz() * b.get_zx();
    c.yy_ = a.get_yx() * b.get_xy() + a.get_yy() * b.get_yy() +
            a.get_yz() * b.get_zy();
    c.yz_ = a.get_yx() * b.get_xz() + a.get_yy() * b.get_yz() +
            a.get_yz() * b.get_zz();
    c.zx_ = a.get_zx() * b.get_xx() + a.get_zy() * b.get_yx() +
            a.get_zz() * b.get_zx();
    c.zy_ = a.get_zx() * b.get_xy() + a.get_zy() * b.get_yy() +
            a.get_zz() * b.get_zy();
    c.zz_ = a.get_zx() * b.get_xz() + a.get_zy() * b.get_yz() +
            a.get_zz() * b.get_zz();
  }

 public:  // creation
  /// @brief Default constructor
  inline Transform()
      : xx_(1.0),
        yx_(0.0),
        zx_(0.0),
        px_(0.0),
        xy_(0.0),
        yy_(1.0),
        zy_(0.0),
        py_(0.0),
        xz_(0.0),
        yz_(0.0),
        zz_(1.0),
        pz_(0.0) {}

  inline Transform(Matrix const& r, Point const& t) {
    rotation(r);
    translation(t);
  }

 public:
  // Accessors

  float get_xx() const { return xx_; }

  float xy() const { return xy_; }

  float xz() const { return xz_; }

  float yx() const { return yx_; }

  float yy() const { return yy_; }

  float yz() const { return yz_; }

  float zx() const { return zx_; }

  float zy() const { return zy_; }

  float zz() const { return zz_; }

  float px() const { return px_; }

  float py() const { return py_; }

  float pz() const { return pz_; }

  float get_xx() { return xx_; }

  float xy() { return xy_; }

  float xz() { return xz_; }

  float yx() { return yx_; }

  float yy() { return yy_; }

  float yz() { return yz_; }

  float zx() { return zx_; }

  float zy() { return zy_; }

  float zz() { return zz_; }

  float px() { return px_; }

  float py() { return py_; }

  float pz() { return pz_; }

  Vector xaxis() const { return Vector(xx_, xy_, xz_); }

  Vector yaxis() const { return Vector(yx_, yy_, yz_); }

  Vector zaxis() const { return Vector(zx_, zy_, zz_); }

  Vector translation() const { return Vector(px_, py_, pz_); }

 public:
  inline Matrix rotation() const {
    return Matrix(xx_, xy_, xz_, yx_, yy_, yz_, zx_, zy_, zz_);
  }

  inline void rotation(Matrix const& m) {
    xx_ = m.get_xx();
    xy_ = m.get_xy();
    xz_ = m.get_xz();
    yx_ = m.get_yx();
    yy_ = m.get_yy();
    yz_ = m.get_yz();
    zx_ = m.get_zx();
    zy_ = m.get_zy();
    zz_ = m.get_zz();
  }

  inline void translation(Vector const& v) {
    px_ = v.get_x();
    py_ = v.get_y();
    pz_ = v.get_z();
  }

 private:
  float xx_, yx_, zx_, px_;
  float xy_, yy_, zy_, py_;
  float xz_, yz_, zz_, pz_;
};

}  // namespace math

#endif /* defined(__REDESIGNC__Transform__) */
