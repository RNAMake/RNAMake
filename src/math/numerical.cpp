//
//  Numeric_Test.cpp
//  REDESIGNC
//
//  Created by Joseph Yesselman on 10/4/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

// RNAMake Headers
#include <base/types.hpp>
#include <math/numerical.hpp>

namespace math {

/// @brief - checks if two float values are equal
int are_floats_equal(double const a, double const b, double tol) {
  if (std::abs(a - b) < tol) {
    return 1;
  } else {
    return 0;
  }
}

/// @brief - checks if two vectors are equal
int are_vector3_equal(Vector3 const &vec, Vector3 const &correct_vec,
                        float tol) {
  if (are_floats_equal(vec.get_x(), correct_vec.get_x(), tol) &&
      are_floats_equal(vec.get_y(), correct_vec.get_y(), tol) &&
      are_floats_equal(vec.get_z(), correct_vec.get_z(), tol)) {
    return 1;
  } else {
    return 0;
  }
}

/// @brief - checks if two points are equal
int are_points_equal(Vector3 const &p1, Vector3 const &p2, float tol) {
  if (are_floats_equal(p1.get_x(), p2.get_x(), tol) &&
      are_floats_equal(p1.get_y(), p2.get_y(), tol) &&
      are_floats_equal(p1.get_z(), p2.get_z(), tol)) {
    return 1;
  } else {
    return 0;
  }
}

/// @brief - checks if two matrices are equal
int are_matrices_equal(Matrix3x3 const &m, Matrix3x3 const &mc) {
  if (!are_floats_equal(m.get_xx(), mc.get_xx()) ||
      !are_floats_equal(m.get_xy(), mc.get_xy()) ||
      !are_floats_equal(m.get_xz(), mc.get_xz()) ||
      !are_floats_equal(m.get_yx(), mc.get_yx()) ||
      !are_floats_equal(m.get_yy(), mc.get_yy()) ||
      !are_floats_equal(m.get_yz(), mc.get_yz()) ||
      !are_floats_equal(m.get_zx(), mc.get_zx()) ||
      !are_floats_equal(m.get_zy(), mc.get_zy()) ||
      !are_floats_equal(m.get_zz(), mc.get_zz())) {
    return 0;
  }

  return 1;
}

/// @brief - checks if two double numbers are within the tolerance of each other
template <>
bool roughly_equal<double>(double const &v1, double const &v2,
                           double tolerance) {
  // TODO maybe add nan checking??
  return std::abs(v1 - v2) < tolerance;
}

/// @brief - checks if two float numbers are within the tolerance of each other
template <>
bool roughly_equal<float>(float const &v1, float const &v2, double tolerance) {
  return std::abs(v1 - v2) < tolerance;
}

/// @brief - checks if two matrices are within the tolerance of each other
template <>
bool roughly_equal<Matrix3x3>(Matrix3x3 const &m1, Matrix3x3 const &m2,
                              double tolerance) {
  return roughly_equal(m1.get_xx(), m2.get_xx(), tolerance) &&
         roughly_equal(m1.get_xy(), m2.get_xy(), tolerance) &&
         roughly_equal(m1.get_xz(), m2.get_xz(), tolerance) &&
         roughly_equal(m1.get_yx(), m2.get_yx(), tolerance) &&
         roughly_equal(m1.get_yy(), m2.get_yy(), tolerance) &&
         roughly_equal(m1.get_yz(), m2.get_yz(), tolerance) &&
         roughly_equal(m1.get_zx(), m2.get_zx(), tolerance) &&
         roughly_equal(m1.get_zy(), m2.get_zy(), tolerance) &&
         roughly_equal(m1.get_zz(), m2.get_zz(), tolerance);
}

/// @brief - checks if two vectors are within the tolerance of each other
template <>
bool roughly_equal<Vector3>(Vector3 const &v1, Vector3 const &v2,
                            double tolerance) {
  return roughly_equal(v1.get_x(), v2.get_x(), tolerance) &&
         roughly_equal(v1.get_y(), v2.get_y(), tolerance) &&
         roughly_equal(v1.get_z(), v2.get_z(), tolerance);
}

} // namespace math