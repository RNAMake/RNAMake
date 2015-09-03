//
//  motif_state_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_unittest__
#define __RNAMake__motif_state_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifStateUnittest : public Unittest {
public:
    
    MotifStateUnittest() {}
    
    ~MotifStateUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_copy();
    
    int
    test_to_str();
    
    int
    test_align();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_state_unittest__) */
