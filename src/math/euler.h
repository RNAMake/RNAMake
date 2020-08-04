//
//  euler.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_euler_h
#define RNAMake_euler_h

#include "math/xyz_vector.h"
#include "math/xyz_matrix.h"

namespace math {

const double _EPS = 2.22044604925e-16 * 4.0;

//assumes 3x3 matrices
void
calc_euler(
        Matrix & M,
        Vector & euler);

struct AxisAngle {
    float angle;
    Point axis;
};

void
axis_angle_from_matrix(
        Matrix &,
        AxisAngle &);

float
degrees(
        float);

}

#endif

