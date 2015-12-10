//
//  residue_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_unittest__
#define __RNAMake__residue_unittest__

#include <stdio.h>
#include <iostream>

//RNAMake Headers
#include "structure/residue.h"
#include "structure/residue_type_set.h"

namespace unittests {

class ResidueUnittest : public Unittest {
public:
    
    int
    test_bead_creation();
    
    int
    test_str_to_residue();
    
    int
    test_get_atom();
    
    int
    test_connected_to();
    
    int
    test_get_beads();
    
    int
    test_copy();
    
    int
    test_to_str();
    
    int
    test_equals();
    
    int
    test_memory_management();
    
public:
    
    int
    run();
    
    int
    run_all();
    
    
};

    
}

#endif /* defined(__RNAMake__residue_unittest__) */
