//
//  basepair.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure/basepair.h"

Basepair::Basepair(
        UuidOP const & res1_uuid,
        UuidOP const & res2_uuid,
        Matrix const & r,
        Point const & d,
        Points const & sugars,
        StringOP const & name,
        StringOP const & x3dna_bp_type,
        primitives::Basepair::BasepairType const & bp_type,
        UuidOP const & uuid) :
        primitives::Basepair(uuid),
        res1_uuid_(res1_uuid),
        res2_uuid_(res2_uuid),
        r_(r),
        d_(d),
        sugars_(sugars),
        name_(name),
        x3dna_bp_type_(x3dna_bp_type),
        bp_type_(bp_type) {}

Basepair::Basepair(
        String const & s,
        UuidOP const & res1_uuid,
        UuidOP const & res2_uuid):
        primitives::Basepair(),
        res1_uuid_(res1_uuid),
        res2_uuid_(res2_uuid) {

    auto spl = split_str_by_delimiter(s, ";");
    d_       = vector_from_str(spl[0]);
    r_       = matrix_from_str(spl[1]);
    sugars_  = vectors_from_str(spl[2]);
    name_    = std::make_shared<String>(spl[3]);
    x3dna_bp_type_ = std::make_shared<String>(spl[4]);
    bp_type_ = static_cast<primitives::Basepair::BasepairType>(std::stoi(spl[5]));
}

Basepair::Basepair(Basepair const & bp):
        primitives::Basepair(bp.uuid_),
        res1_uuid_(bp.res1_uuid_),
        res2_uuid_(bp.res2_uuid_),
        r_(bp.r_),
        d_(bp.d_),
        sugars_(bp.sugars_),
        name_(bp.name_),
        x3dna_bp_type_(bp.x3dna_bp_type_),
        bp_type_(bp.bp_type_) {}

Basepair::Basepair(
        Basepair const & bp,
        UuidOP const & res1_uuid,
        UuidOP const & res2_uuid,
        UuidOP const & bp_uuid):
        primitives::Basepair(bp_uuid),
        res1_uuid_(res1_uuid),
        res2_uuid_(res2_uuid),
        r_(bp.r_),
        d_(bp.d_),
        sugars_(bp.sugars_),
        name_(bp.name_),
        x3dna_bp_type_(bp.x3dna_bp_type_),
        bp_type_(bp.bp_type_)  {

    if(bp_uuid == nullptr) { uuid_ = std::make_shared<Uuid>(); }
}

String
Basepair::to_str() const {
    auto s = vector_to_str(d_) + ";" + matrix_to_str(r_) + ";" + vectors_to_str(sugars_) + ";";
    s += *name_ + ";" + *x3dna_bp_type_ + ";" + std::to_string(bp_type_) + ";";
    return s;
}

state::BasepairOP
Basepair::get_state() {
    auto bp_state = std::make_shared<state::Basepair>(res1_uuid_, res2_uuid_, r_, d_, sugars_, name_,
                                                      x3dna_bp_type_, bp_type_, uuid_);
    return bp_state;
}


bool
are_basepairs_equal(
        BasepairOP const & bp1,
        BasepairOP const & bp2,
        int check_uuids) {

    if(!are_xyzVector_equal(bp1->d(), bp2->d())) { return false; }
    if(!are_xyzMatrix_equal(bp1->r(), bp2->r())) { return false; }
    if(!are_xyzVector_equal(bp1->res1_sugar(), bp2->res1_sugar())) { return false; }
    if(!are_xyzVector_equal(bp1->res2_sugar(), bp2->res2_sugar())) { return false; }
    if((*bp1->name()) != (*bp2->name())) { return false; }
    if((*bp1->x3dna_bp_type()) != (*bp2->x3dna_bp_type())) { return false; }

    if(check_uuids) {
        if((*bp1->uuid()) != (*bp2->uuid())) { return false; }
        if((*bp1->res1_uuid()) != (*bp2->res1_uuid())) { return false; }
        if((*bp1->res2_uuid()) != (*bp2->res2_uuid())) { return false; }
    }

    return true;

}


Point
_calc_center(ResidueOPs const & res) {
    auto center = Point(0, 0, 0);
    auto total = 0.0f;
    for (auto const & r : res) {
        for (auto const & a : *r) {
            if (a == nullptr) { continue; }
            center += a->coords();
            total += 1.0f;
        }
    }
    center /= total;
    return center;
}



