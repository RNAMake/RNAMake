//
//  motif_tree_topology_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_topology_unittest__
#define __RNAMake__motif_tree_topology_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifTreeTopologyUnittest : public Unittest {
public:
    
    MotifTreeTopologyUnittest() {}
    
    ~MotifTreeTopologyUnittest() {}
    
public:
    
    int
    test_creation();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__motif_tree_topology_unittest__) */
