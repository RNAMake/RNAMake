//
//  resource_manager_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "resource_manager_unittest.h"

#include "resources/resource_manager.h"

int
ResourceManagerUnittest::test_get_motif() {
    auto m = ResourceManager::getInstance().get_motif("HELIX.IDEAL");
    if(m->name() != "HELIX.IDEAL") { return 0; }
    return 1;
}


int
ResourceManagerUnittest::run() {
    if (test_get_motif() == 0)        { std::cout << "test_get_motif failed" << std::endl; }
  
    return 0;
}