//
// Created by Joseph Yesselman on 2/14/18.
//

#include <math/euler.h>

namespace math {

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

void
axis_angle_from_matrix(
        Matrix & m,
        AxisAngle & aa) {
    auto epsilon = 0.01; // margin to allow for rounding errors
    auto epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees
    if(std::abs(m.xy() - m.yx()) < epsilon && std::abs(m.xz() - m.zx()) < epsilon && std::abs(m.yz() - m.zy()) < epsilon) {
        // singularity found
        // first check for identity matrix which must have +1 for all terms
        //  in leading diagonaland zero in other terms
        if(std::abs(m.xy() + m.yx()) < epsilon2 && std::abs(m.xz() + m.zx()) < epsilon2 &&
           std::abs(m.yz() + m.zy()) < epsilon2 && std::abs(m.xx() + m.yy() + m.zz() - 3) < epsilon2) {
            // this singularity is identity matrix so angle = 0
            aa.angle = 0;
            aa.axis.x(1); aa.axis.y(0); aa.axis.z(0);
            return;
        }
        aa.angle = M_PI;
        auto xx = (m.xx() + 1) / 2;
        auto yy = (m.yy() + 1) / 2;
        auto zz = (m.zz() + 1) / 2;
        auto xy = (m.xy() + m.yx()) / 4;
        auto xz = (m.xz() + m.zx()) / 4;
        auto yz = (m.yz() + m.zy()) / 4;
        if     ( xx > yy && xx > zz) {
            if(xx < epsilon) { aa.axis.x(0); aa.axis.y(0.7071); aa.axis.z(0.7071); }
            else             { aa.axis.x(sqrt(xx)); aa.axis.y(xy/aa.axis.x()); aa.axis.z(yz/aa.axis.x()); }
        }
        else if(yy > zz ) {
            if(yy < epsilon) { aa.axis.x(0.7071); aa.axis.y(0); aa.axis.z(0.7071);}
            else             { aa.axis.y(sqrt(yy)); aa.axis.x(xy/aa.axis.y()); aa.axis.z(yz/aa.axis.y());}
        }
        else {
            if(zz < epsilon) { aa.axis.x(0.7071); aa.axis.y(0.7071); aa.axis.z(0); }
            else             { aa.axis.z(sqrt(zz)); aa.axis.x(xz/aa.axis.z()); aa.axis.y(yz/aa.axis.z()); }
        }
        return;
    }

    auto s = sqrt((m.zy() - m.yz())*(m.zy() - m.yz()) + (m.xz() - m.zx())*(m.xz() - m.zx()) +
                  (m.yx() - m.xy())*(m.yx() - m.xy()));
    if(std::abs(s) < 0.001) { s = 1; }
    aa.angle = acos((m.xx() + m.yy() + m.zz() - 1) / 2);
    aa.axis.x((m.zy() - m.yz()) / s);
    aa.axis.y((m.xz() - m.xz()) / s);
    aa.axis.z((m.yx() - m.xy()) / s);
}

float
degrees(
        float radians) {
    return radians*(180/M_PI);
}


}