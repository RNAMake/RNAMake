//
// Created by Joseph Yesselman on 12/16/17.
//

#include <iostream>
#include "math/quaternion.h"

namespace math {

std::ostream &
operator<<(
        std::ostream & stream,
        Quaternion const & q) {
    stream << q.get_a() << " " << q.get_b() << " " << q.get_c() << " " << q.get_d();
    return stream;
}

Quaternion
get_random_quaternion() {
    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    auto rand_vals = std::vector<double>(3);
    for (int i = 0; i < 3; i++) { rand_vals[i] = dis(gen); }

    auto r1 = sqrt(1.0 - rand_vals[0]);
    auto r2 = sqrt(rand_vals[0]);
    auto pi2 = M_PI * 2.0;
    auto t1 = pi2 * rand_vals[1];
    auto t2 = pi2 * rand_vals[2];
    return Quaternion(cos(t2) * r2, sin(t1) * r1, cos(t1) * r1, sin(t2) * r2);
}

Matrix
Quaternion::get_rotation_matrix() {
    auto n = dot(*this);
    auto q1 = Quaternion(a_, b_, c_, d_);
    q1 *= sqrt(2.0 / n);

    auto q = std::vector<std::vector<double>>();
    for (int i = 0; i < 4; i++) {
        q.push_back(std::vector<double>(4));
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            q[i][j] = q1[i] * q1[j];
        }
    }

    return Matrix(
            1.0 - q[2][2] - q[3][3], q[1][2] - q[3][0], q[1][3] + q[2][0],
            q[1][2] + q[3][0], 1.0 - q[1][1] - q[3][3], q[2][3] - q[1][0],
            q[1][3] - q[2][0], q[2][3] + q[1][0], 1.0 - q[1][1] - q[2][2]);
}


Quaternion
get_quaternion_from_matrix(
        Matrix const & m) {

    double qx, qy, qz, qw;
    double tr = m.xx() + m.yy() + m.zz();
    if (tr > 0) {
        double s = sqrt(tr + 1.0) * 2;
        qw = 0.25 * s;
        qx = (m.zy() - m.yz()) / s;
        qy = (m.xz() - m.zx()) / s;
        qz = (m.yx() - m.xy()) / s;
    } else if ((m.xx() > m.yy()) && (m.xx() > m.zz())) {
        double s = sqrt(1.0 + m.xx() - m.yy() - m.zz()) * 2;
        qw = (m.zy() - m.yz()) / s;
        qx = 0.25 * s;
        qy = (m.xy() + m.yx()) / s;
        qz = (m.xz() + m.zx()) / s;
    } else if (m.yy() > m.zz()) {
        double s = sqrt(1.0 + m.yy() - m.xx() - m.zz()) * 2;
        qw = (m.xz() - m.zx()) / s;
        qx = (m.xy() + m.yx()) / s;
        qy = 0.25 * s;
        qz = (m.zy() + m.yz()) / s;
    } else {
        double s = sqrt(1.0 + m.zz() - m.xx() - m.yy()) * 2;
        qw = (m.xy() - m.yx()) / s;
        qx = (m.xz() + m.zx()) / s;
        qy = (m.yz() + m.zy()) / s;
        qz = 0.25 * s;
    }

    return Quaternion(qw, qx, qy, qz);
}


void
power_iteration(
        std::vector<std::vector<double>> const & A,
        std::vector<double> & eigen_values,
        int num_simulations) {

    std::random_device rd;
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    auto b_k = std::vector<double>(A.size());
    auto b_k1 = std::vector<double>(A.size());
    double b_k1_norm = 0.0;
    for (int i = 0; i < 2; i++) { b_k[i] = dis(gen); }

    for (int i = 0; i < num_simulations; i++) {
        // calculate the matrix-by-vector product Ab
        dot_vector(A, b_k, b_k1);

        // calculate the norm
        b_k1_norm = norm(b_k1);

        for (int j = 0; j < b_k1.size(); j++) {
            b_k[j] = b_k1[j] / b_k1_norm;
        }
    }

    eigen_values = b_k;

}

}
































