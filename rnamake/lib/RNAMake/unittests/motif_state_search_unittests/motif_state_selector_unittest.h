//
//  motif_state_selector_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/13/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_selector_unittest__
#define __RNAMake__motif_state_selector_unittest__

#include <stdio.h>


#include "unittest.h"

class MotifStateSelectorUnittest : public Unittest {
public:
    MotifStateSelectorUnittest() {}
    
    ~MotifStateSelectorUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_get_children_ms();
    
    int
    test_default_selector();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_state_selector_unittest__) */
