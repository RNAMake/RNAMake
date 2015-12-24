//
//  motif_tree_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMake Headers
#include "motif_data_structures_unittests/motif_tree_unittest.h"

int main(int argc, const char * argv[]) {
    unittests::motif_structures::MotifTreeUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "motif_data_structures_unittests/motif_tree_unittest.h: ";
    std::cout << passed << " PASSED!" << std::endl;
    return 0;
    
}