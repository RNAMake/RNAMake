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
    double cy = sqrt(M.get_zz()*M.get_zz() + M.get_yz()*M.get_yz());
    if(cy > _EPS) {
        //ax = math.atan2( M[k, j],  M[k, k])
        //ay = math.atan2(-M[k, i],  cy)
        //az = math.atan2( M[j, i],  M[i, i])
        euler[0] = atan2( M.get_xy(), M.get_xx());
        euler[1] = atan2(-M.get_xz(), cy);
        euler[2] = atan2( M.get_yz(), M.get_zz());
    }
    else {
        //ax = math.atan2(-M[j, k],  M[j, j])
        //ay = math.atan2(-M[k, i],  cy)
        //az = 0.0
        euler[0] = atan2( M.get_yx(), M.get_yy());
        euler[1] = atan2(-M.get_xz(), cy);
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
    if(std::abs(m.get_xy() - m.get_yx()) < epsilon && std::abs(m.get_xz() - m.get_zx()) < epsilon && std::abs(m.get_yz() - m.get_zy()) < epsilon) {
        // singularity found
        // first check for identity matrix which must have +1 for all terms
        //  in leading diagonaland zero in other terms
        if(std::abs(m.get_xy() + m.get_yx()) < epsilon2 && std::abs(m.get_xz() + m.get_zx()) < epsilon2 &&
           std::abs(m.get_yz() + m.get_zy()) < epsilon2 && std::abs(m.get_xx() + m.get_yy() + m.get_zz() - 3) < epsilon2) {
            // this singularity is identity matrix so angle = 0
            aa.angle = 0;
            aa.axis.set_x(1); aa.axis.set_y(0); aa.axis.set_z(0);
            return;
        }
        aa.angle = M_PI;
        auto xx = (m.get_xx() + 1) / 2;
        auto yy = (m.get_yy() + 1) / 2;
        auto zz = (m.get_zz() + 1) / 2;
        auto xy = (m.get_xy() + m.get_yx()) / 4;
        auto xz = (m.get_xz() + m.get_zx()) / 4;
        auto yz = (m.get_yz() + m.get_zy()) / 4;
        if     ( xx > yy && xx > zz) {
            if(xx < epsilon) { aa.axis.set_x(0); aa.axis.set_y(0.7071); aa.axis.set_z(0.7071); }
            else             { aa.axis.set_x(sqrt(xx)); aa.axis.set_y(xy/aa.axis.get_x()); aa.axis.set_z(yz/aa.axis.get_x()); }
        }
        else if(yy > zz ) {
            if(yy < epsilon) { aa.axis.set_x(0.7071); aa.axis.set_y(0); aa.axis.set_z(0.7071);}
            else             { aa.axis.set_y(sqrt(yy)); aa.axis.set_x(xy/aa.axis.get_y()); aa.axis.set_z(yz/aa.axis.get_y());}
        }
        else {
            if(zz < epsilon) { aa.axis.set_x(0.7071); aa.axis.set_y(0.7071); aa.axis.set_z(0); }
            else             { aa.axis.set_z(sqrt(zz)); aa.axis.set_x(xz/aa.axis.get_z()); aa.axis.set_y(yz/aa.axis.get_z()); }
        }
        return;
    }

    auto s = sqrt((m.get_zy() - m.get_yz())*(m.get_zy() - m.get_yz()) + (m.get_xz() - m.get_zx())*(m.get_xz() - m.get_zx()) +
                  (m.get_yx() - m.get_xy())*(m.get_yx() - m.get_xy()));
    if(std::abs(s) < 0.001) { s = 1; }
    aa.angle = acos((m.get_xx() + m.get_yy() + m.get_zz() - 1) / 2);
    aa.axis.set_x((m.get_zy() - m.get_yz()) / s);
    aa.axis.set_y((m.get_xz() - m.get_xz()) / s);
    aa.axis.set_z((m.get_yx() - m.get_xy()) / s);
}

float
degrees(
        float radians) {
    return radians*(180/M_PI);
}


}