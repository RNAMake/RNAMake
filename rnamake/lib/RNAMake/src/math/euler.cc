//
// Created by Joseph Yesselman on 2/14/18.
//

#include <math/euler.h>

void
calc_euler(
        Matrix & M,
        Vector & euler) {
    // assuming rxyz
    // i = z
    // j = y
    // k = x

    //cy = math.sqrt(M[i, i]*M[i, i] + M[j, i]*M[j, i])
    double cy = sqrt(M.zz()*M.zz() + M.yz()*M.yz());
    if(cy > _EPS) {
        //ax = math.atan2( M[k, j],  M[k, k])
        //ay = math.atan2(-M[k, i],  cy)
        //az = math.atan2( M[j, i],  M[i, i])
        euler[0] = atan2( M.xy(), M.xx());
        euler[1] = atan2(-M.xz(), cy);
        euler[2] = atan2( M.yz(), M.zz());
    }
    else {
        //ax = math.atan2(-M[j, k],  M[j, j])
        //ay = math.atan2(-M[k, i],  cy)
        //az = 0.0
        euler[0] = atan2( M.yx(), M.yy());
        euler[1] = atan2(-M.xz(), cy);
        euler[2] = 0.0;
    }
    std::swap(euler[0], euler[2]);
    euler[0] = -euler[0];
    euler[1] = -euler[1];
    euler[2] = -euler[2];


}
