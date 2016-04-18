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
are_atoms_equal(
    AtomOP const & a1,
    AtomOP const & a2) {
    
    return are_xyzVector_equal(a1->coords(), a2->coords()) &&
           a1->name() == a2->name();
}


