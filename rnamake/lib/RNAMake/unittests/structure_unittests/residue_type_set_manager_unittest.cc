//
//  resource_manager_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/10/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue_type_set_manager_unittest.h"
#include "structure/residue_type_set_manager.h"

int
ResidueTypeSetManagerUnittest::test_creation() {
    //turn on a print statement should only see it created once!
    for(int i = 0; i < 1000; i++) {
        ResidueTypeSetManager::getInstance().residue_type_set();
    }
    
    return 1;
}

int
ResidueTypeSetManagerUnittest::run() {
    
    if (test_creation() == 0)    { std::cout << "test_creation failed" << std::endl; }
    return 0;
    
    
}
