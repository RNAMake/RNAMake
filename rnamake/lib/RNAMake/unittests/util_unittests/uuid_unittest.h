//
//  uuid_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__uuid_unittest__
#define __RNAMake__uuid_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "util/uuid.h"

class UuidUnittest : public Unittest {
public:
    
    int
    test_compare();
    
    int
    test_map();
    
public:
    
    int
    run();
    
    
};

#endif /* defined(__RNAMake__uuid_unittest__) */
