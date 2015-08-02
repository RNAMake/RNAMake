//
//  resource_manager_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/10/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__resource_manager_unittest__
#define __RNAMake__resource_manager_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class ResidueTypeSetManagerUnittest : public Unittest {
public:
    
    ResidueTypeSetManagerUnittest() {}
    
    ~ResidueTypeSetManagerUnittest() {}
    
public:
    
    int
    test_creation();
    
public:
    
    int
    run();
    
private:
    
    
};

#endif /* defined(__RNAMake__resource_manager_unittest__) */
