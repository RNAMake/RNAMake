//
//  tree_memory_test.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMAke Headers
#include "data_structure/tree/tree.h"
#include "util/random_number_generator.h"

int main(int argc, const char * argv[]) {
    RandomNumberGenerator rng;
    
    for(int i = 0; i < 100000; i++) {
        TreeDynamic<int> t;
        int max = rng.randrange(100);
        for(int j = 0; j < max; j++) {
            t.add_data(rng.randrange(100));
        }
    }

}