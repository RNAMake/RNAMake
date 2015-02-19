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
    designable_ = std::map<Uuid, int> ();
}

Pose::Pose(
    StructureOP const & structure,
    BasepairOPs const & basepairs) {
    structure_ = structure;
    basepairs_ = basepairs;
    setup_basepair_ends();
    designable_ = std::map<Uuid, int> ();
}



void
Pose::set_bp_designable(BasepairOP const & bp) {
    designable_[bp->uuid()] = 1;
}