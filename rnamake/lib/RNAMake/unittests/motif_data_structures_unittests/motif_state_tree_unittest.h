//
//  motif_state_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_tree_unittest__
#define __RNAMake__motif_state_tree_unittest__

#include <stdio.h>


//RNAMake Headers
#include "unittest.h"

class MotifStateTreeUnittest : public Unittest {
public:
    
    MotifStateTreeUnittest() {}
    
    ~MotifStateTreeUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_add_state();
    
    int
    test_from_mt();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_state_tree_unittest__) */
