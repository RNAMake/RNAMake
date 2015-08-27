//
//  pose_factory_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pose_factory_unittest__
#define __RNAMake__pose_factory_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class PoseFactoryUnittest : public Unittest {
public:
    
    PoseFactoryUnittest() {}
    
    ~PoseFactoryUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_load_pose();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__pose_factory_unittest__) */
