//
// Created by Joseph Yesselman on 1/29/17.

#include "base/util.hpp"
#include "math/numerical.h"
#include "motif_state/basepair.h"

namespace state {

Basepair::Basepair(
        String const & s,
        Uuid const & res1_uuid,
        Uuid const & res2_uuid,
        Uuid const & uuid):
        primitives::Basepair(res1_uuid, res2_uuid, uuid) {

    auto spl = split_str_by_delimiter(s, ";");
    d_ = vector_from_str(spl[0]);
    r_ = matrix_from_str(spl[1]);
    auto sugars = vectors_from_str(spl[2]);
    res1_sugar_ = sugars[0];
    res2_sugar_ = sugars[1];
    x3dna_bp_type_ = static_cast<X3dna::X3dnaBPType>(std::stoi(spl[4]));
    bp_type_ = static_cast<primitives::Basepair::BasepairType>(std::stoi(spl[5]));
    name_ = std::make_shared<SimpleString>(spl[3]);

}

String
Basepair::to_str() {
    auto sugars = Points{res1_sugar_, res2_sugar_};
    auto s = vector_to_str(d_) + ";" + matrix_to_str(r_) + ";" + vectors_to_str(sugars) + ";";
    s     += name_->to_str() + ";" + std::to_string(x3dna_bp_type_) + ";";
    s     += std::to_string(bp_type_) + ";";
    return s;
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
    if(bp1->name() != bp2->name()) { return false; }
    if(bp1->x3dna_bp_type() != bp2->x3dna_bp_type()) { return false; }

    if(check_uuids) {
        if(bp1->uuid() != bp2->uuid()) { return false; }
        if(bp1->res1_uuid() != bp2->res1_uuid()) { return false; }
        if(bp1->res2_uuid() != bp2->res2_uuid()) { return false; }
    }

    return true;

}

}