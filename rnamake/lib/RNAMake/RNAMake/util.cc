//
//  util.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util.h"
#include "types.h"
#include "atom.h"

Point
center(
    AtomOPs const & atoms) {
    Point center(0, 0, 0);
    for(auto const & a : atoms) {
        if(a == NULL) {
            continue;
        }
        center += a->coords();
    }
    
    return center / float(atoms.size());
}
