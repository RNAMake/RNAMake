//
//  graph_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__graph_unittest__
#define __RNAMake__graph_unittest__

#include <stdio.h>

#include "unittest.h"

class GraphUnittest : public Unittest {
public:
    GraphUnittest() {}
    
    ~GraphUnittest() {}
    
public:
    
    int
    test_nodes();
    
    int
    test_creation();
    
    int
    test_add();
    
    int
    test_connect();
    
    int
    test_remove();
    
    int
    test_iteration();
    
public:
    
    int
    run();
    
    void
    run_all();
    
};


#endif /* defined(__RNAMake__graph_unittest__) */
