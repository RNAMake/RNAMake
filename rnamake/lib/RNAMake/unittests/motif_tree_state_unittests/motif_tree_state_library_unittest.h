//
//  motif_tree_state_library_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/15/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_library_unittest__
#define __RNAMake__motif_tree_state_library_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifTreeStateLibraryUnittest : public Unittest {
public:
    
    MotifTreeStateLibraryUnittest() {}
    
    ~MotifTreeStateLibraryUnittest() {}
    
public:
    
    int
    test_creation();
    
public:
    
    int
    run();
    
    void
    run_all();
    
};



#endif /* defined(__RNAMake__motif_tree_state_library_unittest__) */
