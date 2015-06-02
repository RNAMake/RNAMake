//
//  tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/26/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__tree_unittest__
#define __RNAMake__tree_unittest__

#include <stdio.h>

#include "unittest.h"

class TreeUnittest : public Unittest {
public:
    TreeUnittest() {}
    
    ~TreeUnittest() {}
    
public:
    
    int
    test_node();
    
    int
    test_creation();
    
    int
    test_remove_node();
    
    int
    test_get_index();
    
    int
    test_iter();
    
public:
    
    int
    run();

    
};

#endif /* defined(__RNAMake__tree_unittest__) */
