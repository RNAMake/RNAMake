//
//  vienna_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>


//RNAMake Headers
#include "vienna_unittests/vienna_unittest.h"

int main(int argc, const char * argv[]) {
    unittests::vienna_unittests::ViennaUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "vienna_unittests/vienna_unittest.h: ";
    std::cout << passed << " PASSED!" << std::endl;
    return 0;
    
}