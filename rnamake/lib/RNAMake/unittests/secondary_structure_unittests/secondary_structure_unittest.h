//
//  secondary_structure_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_unittest__
#define __RNAMake__secondary_structure_unittest__

#include <stdio.h>

#include "unittest.h"

class SecondaryStructureUnittest : public Unittest {
public:
    SecondaryStructureUnittest() {}
    
    ~SecondaryStructureUnittest() {}
    
public:
    
    int
    test_creation();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__secondary_structure_unittest__) */
