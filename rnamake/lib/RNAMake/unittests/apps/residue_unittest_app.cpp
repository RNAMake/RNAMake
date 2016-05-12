//
//  residue_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 5/7/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

//RNAMake Headers
#include "structure_unittests/residue_unittest.h"

int main(int argc, const char * argv[]) {
    unittests::ResidueUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "structure_unittests/residue_unittest.h: ";
    std::cout << passed << " PASSED!" << std::endl;
    return 0;
    
}