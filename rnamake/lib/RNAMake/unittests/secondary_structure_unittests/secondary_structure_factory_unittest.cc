//
//  secondary_structure_factory_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_factory_unittest.h"

#include "secondary_structure/secondary_structure_factory.h"

int
SecondaryStructureFactoryUnittest::test_creation() {
    sstruct::SecondaryStructureFactory sf;
    auto ss = sf.get_structure("GGGGG+CCCCC", "(((((+)))))");
    
    return 0;
}


int
SecondaryStructureFactoryUnittest::run() {
    test_creation();
    return 1;
    
}