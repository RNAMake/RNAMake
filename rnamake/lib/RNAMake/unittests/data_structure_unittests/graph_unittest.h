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
    test_creation();
    
    
public:
    
    int
    run();
    
    
};


#endif /* defined(__RNAMake__graph_unittest__) */
