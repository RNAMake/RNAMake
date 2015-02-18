//
//  pose.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pose__
#define __RNAMake__pose__

#include <stdio.h>
#include <memory>
#include "motif.h"

class Pose : public Motif {
public:
    Pose() {} 
    Pose(MotifOP const &);
};

typedef std::shared_ptr<Pose> PoseOP;

#endif /* defined(__RNAMake__pose__) */

