//
// Created by Joe Yesselman on 5/28/22.
//

#ifndef RNAMAKE_SRC_MATH_ROTATION_HPP_
#define RNAMAKE_SRC_MATH_ROTATION_HPP_

#include <math.h>
#include <math/matrix_3x3.hpp>

namespace math {

/// @brief - multiplies two matrices and spits out the value
inline Matrix3x3 rotation_between_frames(const Matrix3x3 &m1,
                                         const Matrix3x3 &m2) {
  Matrix3x3 rot = m1.get_transposed() * m2;
  rot.unitarize();
  return rot;
}

/// @brief - multiplies two matrices and saves the value in "rot"
inline void rotation_between_frames(const Matrix3x3 &m1, const Matrix3x3 &m2,
                                    Matrix3x3 &rot /* return */) {
  rot = m1.get_transposed() * m2;
  rot.unitarize();
}

/// http://www.boris-belousov.net/2016/12/01/quat-dist/
inline Real difference_between_frames(const Matrix3x3 &p, const Matrix3x3 &q) {
  Matrix3x3 r = p * q.get_transposed();
  Real tr = r.get_trace();
  return acos((tr - 1)/2);
}


} // namespace math

#endif // RNAMAKE_SRC_MATH_ROTATION_HPP_
