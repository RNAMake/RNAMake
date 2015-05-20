//
//  ss_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_tree_unittest__
#define __RNAMake__ss_tree_unittest__

#include <stdio.h>

#include "unittest.h"

class SS_TreeUnittest : public Unittest {
public:
    SS_TreeUnittest() {}
    
    ~SS_TreeUnittest() {}
    
public:
    
    int
    test_creation();

public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__ss_tree_unittest__) */
