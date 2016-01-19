//
//  beads_vs_points.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__beads_vs_points__
#define __RNAMake__beads_vs_points__

#include <stdio.h>
#include "structure/residue.h"
#include "util/random_number_generator.h"

void
inline
point_sterics(
    Points const & points,
    RandomNumberGenerator & rng) {
    
    Point p;
    for(int i = 0; i < 1000000000; i++) {
        int count = 0;
        p = points[rng.randrange((int)points.size())];
        for (auto const & p1 : points) {
            if (p.distance(p1) < 0.5) {
                count += 1;
            }
        }
    }
    
    
}

void
inline
beads_sterics(
    Beads const & beads,
    RandomNumberGenerator & rng) {
    
    Bead b;
    for(int i = 0; i < 1000000000; i++) {
        int count = 0;
        b = beads[rng.randrange((int)beads.size())];
        for (auto const & b1 : beads) {
            if (b.center().distance(b1.center()) < 0.5) {
                count += 1;
            }
        }
    }
    
    
}


#endif /* defined(__RNAMake__beads_vs_points__) */
