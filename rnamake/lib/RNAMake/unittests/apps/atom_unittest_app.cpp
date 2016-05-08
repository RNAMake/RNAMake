//
//  atom_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/7/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMake Headers
#include "structure_unittests/atom_unittest.h"

int main(int argc, const char * argv[]) {
    unittests::AtomUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "structure_unittests/atom_unittest.h: ";
    std::cout << passed << " PASSED!" << std::endl;
    return 0;
    
}