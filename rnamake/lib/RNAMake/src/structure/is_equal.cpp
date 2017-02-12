//
//  is_equal.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "is_equal.hpp"
#include "math/numerical.h"







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