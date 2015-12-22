//
//  cl_option_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/20/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMake Headers
#include "base_unittests/cl_option_unittest.h"

int main(int argc, const char * argv[]) {
    CL_OptionUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "base_unittests/cl_option_unittest: " << passed << " PASSED!" << std::endl;
    return 0;
    
}