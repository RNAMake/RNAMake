//
// Created by Joseph Yesselman on 1/29/17.
//

#ifndef RNAMAKE_MOTIF_STATE_BASEPAIR_H
#define RNAMAKE_MOTIF_STATE_BASEPAIR_H

#include "util/uuid.h"
#include "primitives/basepair.h"
#include "math/xyz_matrix.h"
#include "math/transform.h"
#include "util/x3dna.h"

namespace state {

struct TransformInfo {
    TransformInfo():
        r(Matrix()),
        d(Point()),
        sugars(Points(2)) {}

    Matrix r;
    Point d;
    Points sugars;
};

class Basepair : public primitives::Basepair {
public:
    inline
    Basepair(
            Uuid const & res1_uuid,
            Uuid const & res2_uuid,
            Matrix const & r,
            Point const & d,
            Point const & res1_sugar,
            Point const & res2_sugar,
            Chars const & name,
            X3dna::X3dnaBPType const x3dna_bp_type,
            primitives::Basepair::BasepairType const bp_type,
            Uuid const & uuid) :
            primitives::Basepair(uuid),
            res1_uuid_(res1_uuid),
            res2_uuid_(res2_uuid),
            r_(r),
            d_(d),
            res1_sugar_(res1_sugar),
            res2_sugar_(res2_sugar),
            name_(name),
            x3dna_bp_type_(x3dna_bp_type),
            bp_type_(bp_type) {}

    inline
    Basepair(Basepair const & bp):
            primitives::Basepair(bp.uuid_),
            res1_uuid_(bp.res1_uuid_),
            res2_uuid_(bp.res2_uuid_),
            r_(bp.r_),
            d_(bp.d_),
            res1_sugar_(bp.res1_sugar_),
            res2_sugar_(bp.res2_sugar_),
            name_(bp.name_),
            x3dna_bp_type_(bp.x3dna_bp_type_),
            bp_type_(bp.bp_type_) {}

    Basepair(
            String const &,
            Uuid const &,
            Uuid const &,
            Uuid const &);

    ~Basepair() {}

public:
    inline
    void
    move(Point const & p) {
        res1_sugar_ += p;
        res2_sugar_ += p;
        d_ += p;
    }

    inline
    void
    transform(Transform const & t) {
        auto r_T = t.rotation().transpose();
        auto new_r = Matrix();
        dot(r_, r_T, new_r);
        r_ = new_r;
        res1_sugar_ = dot_vector(r_T, res1_sugar_) + t.translation();
        res2_sugar_ = dot_vector(r_T, res2_sugar_) + t.translation();
        d_ = dot_vector(r_T, d_) + t.translation();
    }

    String
    to_str();

public:

    void
    get_transforming_r_and_t (
            Basepair const & o_state, //state with desired rotation and translation
            TransformInfo & r_state) {

        auto r_T = Matrix();
        transpose(r_, r_T);
        auto diff_sugars = Points(2);

        //calculate transforming rotation matrix and store it in r_state (result state)
        dot(r_T, o_state.r_, r_state.r);
        r_state.r.unitarize();

        //rotate sugars to new position and store them in r_state
        transpose(r_state.r, r_T);

        dot_vector(r_T, o_state.res1_sugar_, r_state.sugars[0]);
        dot_vector(r_T, o_state.res2_sugar_, r_state.sugars[1]);

        auto diff = -o_state.d_ + d_;
        int i;

        for(i = 0; i < 2; i++) {
            r_state.sugars[i] += diff;
        }

        diff_sugars[0] = res1_sugar_ - r_state.sugars[0];
        diff_sugars[1] = res2_sugar_ - r_state.sugars[1];

        diff = (diff_sugars[0] + diff_sugars[1]) / 2.0f;
        r_state.d = -o_state.d_ + diff + d_;
        r_state.r.transpose();
    }

    void
    fast_transform(TransformInfo const & rt) {

        r_ = dot(r_, rt.r);
        r_.unitarize();
        d_ = dot_vector(rt.r, d_) + rt.d;
        res1_sugar_ = dot_vector(rt.r, res1_sugar_) + rt.d;
        res2_sugar_ = dot_vector(rt.r, res2_sugar_) + rt.d;
    }

public: // getters

    inline
    Matrix const &
    r() const { return r_; }

    inline
    Point const &
    d() const { return d_; }

    inline
    Points
    sugars() const { return Points{res1_sugar_, res2_sugar_}; }

    inline
    Point const &
    res1_sugar() const { return res1_sugar_; }

    inline
    Point const &
    res2_sugar() const { return res2_sugar_; }

    inline
    X3dna::X3dnaBPType const
    x3dna_bp_type() const { return x3dna_bp_type_; }

    inline
    Uuid const &
    uuid() const { return uuid_; }

    inline
    Uuid const &
    res1_uuid() const { return res1_uuid_; }

    inline
    Uuid const &
    res2_uuid() const { return res2_uuid_; }

    inline
    Chars const &
    name() const { return name_; }

    inline
    primitives::Basepair::BasepairType const
    bp_type() const { return bp_type_; }


private:
    Matrix r_;
    Point d_, res1_sugar_, res2_sugar_;
    Uuid res1_uuid_, res2_uuid_;
    Chars name_;
    X3dna::X3dnaBPType x3dna_bp_type_;
    primitives::Basepair::BasepairType bp_type_;
};

typedef std::shared_ptr<Basepair> BasepairOP;
typedef std::vector<BasepairOP>   BasepairOPs;

bool
are_basepairs_equal(
        BasepairOP const & bp1,
        BasepairOP const & bp2,
        int check_uuids = 1);

}

#endif //TEST_BASEPAIR_H
