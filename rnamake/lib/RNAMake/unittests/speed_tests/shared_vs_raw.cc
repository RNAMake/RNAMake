//
//  shared_vs_raw.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/random_number_generator.h"
#include "shared_vs_raw.h"

void raw_test(Tree<DerivedNode> & t) {
    std::vector<int> data(1000);
    for(int i = 0; i < 100000; i++) {
        for(int j = 0; j < 1000; j++) {
            data[j] = t.node(j)->i_ + t.node(j)->i_;
        }
    }
}

void shared_test(Tree<DerivedNode> & t) {
    std::vector<int> data(1000);
    for(int i = 0; i < 100000; i++) {
        for(int j = 0; j < 1000; j++) {
            data[j] = t.node_unique(j)->i_ + t.node_unique(j)->i_;
        }
    }
}

int main(int argc, const char * argv[]) {
    Tree<DerivedNode> t;
    RandomNumberGenerator rng;
    for(int i = 0; i < 1000; i++) {
        t.add_data(rng.randrange(1000));
    }
    
    raw_test(t);
    shared_test(t);
    return 0;
}