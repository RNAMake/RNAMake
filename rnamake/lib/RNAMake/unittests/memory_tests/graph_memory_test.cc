//
//  graph_memory_test.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMAke Headers
#include "data_structure/graph/graph.h"
#include "util/random_number_generator.h"

int main(int argc, const char * argv[]) {
    RandomNumberGenerator rng;
    
    for(int i = 0; i < 10000; i++) {
        GraphDynamic<int> g;
        int max = rng.randrange(100);
        for(int j = 0; j < max; j++) {
            g.add_data(rng.randrange(100));
        }
    }
    
    for(int i = 0; i < 10000; i++) {
        GraphStatic<int> g1;
        int max = rng.randrange(100);
        for(int j = 0; j < max; j++) {
            g1.add_data(j, -1, -1, -1, 2);
        }
    }
    
}