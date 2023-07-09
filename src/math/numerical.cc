//
//  Numeric_Test.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/4/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

// RNAMake Headers
#include "math/numerical.h"

#include "base/types.h"
#include "math/xyz_matrix.h"
#include "math/xyz_vector.h"

namespace math {

int are_floats_equal(double const a, double const b, double tol) {
  if (std::abs(a - b) < tol) {
    return 1;
  } else {
    return 0;
  }
}

int are_xyzVector_equal(Vector const& vec, Vector const& correct_vec,
                        float tol) {
  if (are_floats_equal(vec.x(), correct_vec.x(), tol) &&
      are_floats_equal(vec.y(), correct_vec.y(), tol) &&
      are_floats_equal(vec.z(), correct_vec.z(), tol)) {
    return 1;
  } else {
    return 0;
  }
}

int are_xyzVectors_equal(Vectors const& v, Vectors const& vc) {
  if (v.size() != vc.size()) {
    return 0;
  }

  for (int i = 0; i < v.size(); i++) {
    if (!are_xyzVector_equal(v[i], vc[i])) {
      return 0;
    }
  }

  return 1;
}

int are_xyzMatrix_equal(Matrix const& m, Matrix const& mc) {
  return are_floats_equal(m.xx(), mc.xx()) &&
         are_floats_equal(m.xy(), mc.xy()) &&
         are_floats_equal(m.xz(), mc.xz()) &&
         are_floats_equal(m.yx(), mc.yx()) &&
         are_floats_equal(m.yy(), mc.yy()) &&
         are_floats_equal(m.yz(), mc.yz()) &&
         are_floats_equal(m.zx(), mc.zx()) &&
         are_floats_equal(m.zy(), mc.zy()) && are_floats_equal(m.zz(), mc.zz());
}

template <>
bool roughly_equal<double>(double const& v1, double const& v2,
                           double tolerance) {
  // TODO maybe add nan checking??
  return std::abs(v1 - v2) < tolerance;
}

template <>
bool roughly_equal<float>(float const& v1, float const& v2, double tolerance) {
  return std::abs(v1 - v2) < tolerance;
}

template <>
bool roughly_equal<Matrix>(Matrix const& m1, Matrix const& m2,
                           double tolerance) {
  return roughly_equal(m1.xx(), m2.xx(), tolerance) &&
         roughly_equal(m1.xy(), m2.xy(), tolerance) &&
         roughly_equal(m1.xz(), m2.xz(), tolerance) &&
         roughly_equal(m1.yx(), m2.yx(), tolerance) &&
         roughly_equal(m1.yy(), m2.yy(), tolerance) &&
         roughly_equal(m1.yz(), m2.yz(), tolerance) &&
         roughly_equal(m1.zx(), m2.zx(), tolerance) &&
         roughly_equal(m1.zy(), m2.zy(), tolerance) &&
         roughly_equal(m1.zz(), m2.zz(), tolerance);
}

template <>
bool roughly_equal<Vector>(Vector const& v1, Vector const& v2,
                           double tolerance) {
  return roughly_equal(v1.x(), v2.x(), tolerance) &&
         roughly_equal(v1.y(), v2.y(), tolerance) &&
         roughly_equal(v1.z(), v2.z(), tolerance);
}
}  // namespace math
