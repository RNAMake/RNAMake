//
//  beads_vs_points.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "beads_vs_points.h"


int main(int argc, const char * argv[]) {
    Points points(100);
    RandomNumberGenerator rng;
    for(int i = 0; i < 100; i++) {
        points[i].x(rng.rand());
        points[i].y(rng.rand());
        points[i].z(rng.rand());

    }
    
    Beads beads(100);
    int i = 0;
    for(auto const & p : points) {
        beads[i] = Bead(p, PHOS);
        i++;
    }
    
    //point_sterics(points, rng);
    beads_sterics(beads, rng);
    
    
    return 1;

}
