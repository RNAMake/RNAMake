//
//  Numeric_Test.h
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/4/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#ifndef __REDESIGNC__Numeric_Test__
#define __REDESIGNC__Numeric_Test__

#include <iostream>

// RNAMake Headers
#include <math/matrix_3x3.hpp>
#include <math/vector_3.hpp>

namespace math {

int are_floats_equal(double const a, double const b, double tol = 0.001);

int are_points_equal(Vector3 const &p1, Vector3 const &p2, float tol = 0.001);

int are_matrices_equal(Matrix3x3 const &, Matrix3x3 const &);

template <typename T>
bool roughly_equal(T const &v1, T const &v2, double tolerance = 0.001) {
  if (v1.size() != v2.size()) {
    return false;
  }

  const auto it_len = v1.size();

  for (auto ii = 0; ii < it_len; ++ii) {
    if (!roughly_equal(v1[ii], v2[ii], tolerance)) {
      return false;
    }
  }
  return true;
}
// template specialiation for doubles
template <> bool roughly_equal<double>(double const &, double const &, double);

template <> bool roughly_equal<float>(float const &, float const &, double);

template <>
bool roughly_equal<Matrix3x3>(Matrix3x3 const &, Matrix3x3 const &, double);

template <>
bool roughly_equal<Vector3>(Vector3 const &, Vector3 const &, double);

Matrix3x3 get_random_rotation_matrix();

} // namespace math

#endif /* defined(__REDESIGNC__Numeric_Test__) */