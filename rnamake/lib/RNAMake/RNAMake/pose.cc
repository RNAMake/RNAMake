//
//  pose.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "pose.h"

Pose::Pose(MotifOP const & m) {
    structure_ = m->structure();
    basepairs_ = m->basepairs();
    ends_ = m->ends();
    _cache_basepair_frames();
}
