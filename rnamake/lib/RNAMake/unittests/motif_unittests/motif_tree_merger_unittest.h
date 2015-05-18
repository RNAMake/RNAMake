//
//  motif_tree_merger_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_merger_unittest__
#define __RNAMake__motif_tree_merger_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifTreeMergerUnittest : public Unittest {
public:
    
    MotifTreeMergerUnittest() {}
    
    ~MotifTreeMergerUnittest() {}
    
public:
    
    int
    test_merger();
    
public:
    
    int
    run();
    
    void
    run_all();
    
    
};

#endif /* defined(__RNAMake__motif_tree_merger_unittest__) */
