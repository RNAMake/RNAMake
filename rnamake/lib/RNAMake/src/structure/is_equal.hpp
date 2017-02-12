//
//  is_equal.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef is_equal_hpp
#define is_equal_hpp

#include <stdio.h>

#include "structure/structure.h"
#include "structure/basepair.h"



bool
are_structures_equal(
        StructureOP const & s1,
        StructureOP const & s2,
        int check_uuids = 1);

bool
are_basepairs_equal(
        BasepairOP const & bp1,
        BasepairOP const & bp2,
        int check_uuids = 1);

#endif /* is_equal_hpp */
