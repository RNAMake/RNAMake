//
//  structure_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure_unittest__
#define __RNAMake__structure_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "structure/structure.h"


class StructureUnittest : public Unittest {
public:
    
    StructureUnittest();
    
    ~StructureUnittest() {}
    
public:
    
    int
    test_build_chains();
    
    int
    test_move();
    
    int
    test_transform();
    
    int
    test_get_residue();
    
    int
    test_creation_from_pdb();
    
public:
    
    int
    run();
    
    void
    run_all();
    
private:
    Structure s_;
    
    
    
};

#endif /* defined(__RNAMake__structure_unittest__) */
