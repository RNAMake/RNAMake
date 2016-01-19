//
//  build_secondary_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__build_secondary_structure__
#define __RNAMake__build_secondary_structure__

#include <stdio.h>

#include "secondary_structure/pose.h"

namespace unittests {
    
class BuildSecondaryStructure {
public:
    BuildSecondaryStructure() {}
    
    ~BuildSecondaryStructure() {}
    
public:

    sstruct::PoseOP
    build_helix(int size = 10,
                int opt = 0);
    
    sstruct::PoseOP
    build_hairpin(int size = 10,
                  int opt = 0);
    
};
    
}

#endif /* defined(__RNAMake__build_secondary_structure__) */
