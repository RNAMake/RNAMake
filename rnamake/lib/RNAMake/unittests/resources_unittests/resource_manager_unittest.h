//
//  resource_manager_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager_unittest__
#define __RNAMake__resource_manager_unittest__

#include <stdio.h>

#include "unittest.h"

class ResourceManagerUnittest : public Unittest {
public:
    
    ResourceManagerUnittest() {}
    
    ~ResourceManagerUnittest() {}
    
public:
    
    int
    test_get_motif();
    
    
public:
    
    int
    run();
        
};

#endif /* defined(__RNAMake__resource_manager_unittest__) */
