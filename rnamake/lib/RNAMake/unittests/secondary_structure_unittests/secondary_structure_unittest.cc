//
//  secondary_structure_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_unittest.h"

#include "uuid/uuid.h"
#include "secondary_structure/residue.h"
#include "secondary_structure/secondary_structure.h"


int
SecondaryStructureUnittest::test_creation() {

    sstruct::Residue r("G", ".", 1, "A", Uuid());
    String s = r.to_str();
    auto r2 = sstruct::str_to_residue(s);
    
    if(r.num() != r2.num() || r.chain_id() != r2.chain_id()) {
        return 0;
    }
    
    auto r3 = r.copy();
    
    if(!(r.uuid() == r3.uuid())) {
        return 0;
    }
    
    
    
    return 1;
}


int
SecondaryStructureUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    return 1;
}