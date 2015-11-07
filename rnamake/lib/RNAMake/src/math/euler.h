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

const double _EPS = 2.22044604925e-16 * 4.0;

//assumes 3x3 matrices
inline
void
calc_euler(
    Matrix & M,
    Vector & euler) {
    
    double cy = sqrt(M.xx()*M.xx() + M.yx()*M.yx());
    if(cy > _EPS) {
        euler[0] = atan2( M.zy(), M.zz());
        euler[1] = atan2(-M.zx(), cy);
        euler[2] = atan2( M.yx(), M.xx());
    }
    else {
        euler[0] = atan2( M.yz(), M.yy());
        euler[1] = atan2(-M.zx(), cy);
        euler[2] = 0.0;
    }
    for(int i = 0; i < 3; i++){
        if(euler[i] > 6.14) {
            euler[i] -= 6.14;
        }
        if(euler[i] < 0) {
            euler[i] += 6.14;
        }
    }
    
    //'sxyz': (0, 0, 0, 0)
    //_NEXT_AXIS = [1, 2, 0, 1]
    
}


#endif
