//
//  motif_state_search_unittest_app.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>

//RNAMake Headers
#include "motif_state_search_unittests/motif_state_search_unittest.h"

int main(int argc, const char * argv[]) {
    unittests::motif_state_search::MotifStateSearchUnittest test;
    int failed = test.run_all();
    int passed = test.size() - failed;
    std::cout << "motif_state_search_unittests/motif_state_search_unittest.h: ";
    std::cout << passed << " PASSED!" << std::endl;
    return 0;
    
}