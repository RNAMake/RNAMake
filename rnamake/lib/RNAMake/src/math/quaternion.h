//
// Created by Joseph Yesselman on 12/16/17.
//

#ifndef RNAMAKE_NEW_QUATERNION_H
#define RNAMAKE_NEW_QUATERNION_H

#define _USE_MATH_DEFINES

#include <math.h>
#include <random>

#include "math/xyz_matrix.h"

class Quaternion {
public:
    Quaternion():
            a_(0),
            b_(0),
            c_(0),
            d_(0) {}

    Quaternion(
            double val):
            a_(val),
            b_(val),
            c_(val),
            d_(val) {}

    Quaternion(
            double a,
            double b,
            double c,
            double d):
            a_(a),
            b_(b),
            c_(c),
            d_(d) {}

public:
    friend
    std::ostream &
    operator<< (
            std::ostream & ,
            Quaternion const &);

    inline
    double
    operator[] (
            int i) {
        if     (i == 0) { return a_; }
        else if(i == 1) { return b_; }
        else if(i == 2) { return c_; }
        else if(i == 3) { return d_; }
        else            {
            throw std::runtime_error("invalid index");
        }
    }

    inline
    void
    operator += (
            double const v) {
        a_ += v;
        b_ += v;
        c_ += v;
        d_ += v;
    }

    inline
    void
    operator *= (
            double const v) {
        a_ *= v;
        b_ *= v;
        c_ *= v;
        d_ *= v;
    }
public:

public:
    inline
    double
    dot(
            Quaternion const & q) const {
        return a_*q.a_ + b_*q.b_ + c_*q.c_ + d_*q.d_;
    }

    Matrix
    get_rotation_matrix();

public:
    inline
    double
    get_a() const{ return a_; }

    inline
    double
    get_b() const { return b_; }

    inline
    double
    get_c() const { return c_; }

    inline
    double
    get_d() const { return d_; }

private:
    double a_, b_, c_, d_;

};

Quaternion
get_random_quaternion();


inline
void
dot_vector(
        std::vector<std::vector<double>> const & m,
        std::vector<double> const & v,
        std::vector<double> & vr) {

    for(int i = 0; i < m.size(); i++) {
        vr[i] = 0;
        for(int j = 0; j < m[i].size(); j++) {
            vr[i] += m[i][j]*v[j];
        }
    }

}

inline
double
norm(
        std::vector<double> const & v) {
    double norm = 0;
    for(auto const & e : v) {
        norm += e*e;
    }
    return sqrt(norm);
}

void
power_iteration(
        std::vector<std::vector<double>> const &,
        std::vector<double> &,
        int);


#endif //RNAMAKE_NEW_QUATERNION_H


























