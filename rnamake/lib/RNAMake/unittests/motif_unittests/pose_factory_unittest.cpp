//
//  pose_factory_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 8/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/pose_factory.h"
#include "pose_factory_unittest.h"

int
PoseFactoryUnittest::test_creation() {
    PoseFactory pf;
    return 1;
}

int
PoseFactoryUnittest::run() {
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    return 1;
}