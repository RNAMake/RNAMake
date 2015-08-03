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
    if(ss->motifs("BP_STEP").size() != 4 || ss->motifs("ALL").size() != 4) {
        return 0;
    }

    auto ss1 = sf.get_structure("GGGAAGG+CCCCC", "(((..((+)))))");
    if(ss1->motifs("BULGE").size() != 1) {
        return 0;
    }
    
    auto ss2 = ss1->copy();
    if(ss2.motifs("BP_STEP").size() != 3 || ss2.motifs("ALL").size() != 4 ||
       ss2.motifs("BULGE").size() != 1 ) {
        return 0;
    }
    
    return 1;
}


int
SecondaryStructureFactoryUnittest::run() {
    if (test_creation() == 0)          {  std::cout << "test_creation failed" << std::endl; }
    return 1;
    
}