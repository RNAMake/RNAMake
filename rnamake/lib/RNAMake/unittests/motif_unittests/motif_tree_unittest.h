//
//  motif_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_unittest__
#define __RNAMake__motif_tree_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifTreeUnittest : public Unittest {
public:
  
    MotifTreeUnittest() {}
    
    ~MotifTreeUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_add_motif();
    
    int
    test_motif_tree_to_str();
    
    int
    test_remove_node();
    
    int
    test_remove_node_level();
    
public:
    
    int
    run();
    
    void
    run_all();
    
    
};

#endif /* defined(__RNAMake__motif_tree_unittest__) */
