//
//  pose_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pose_factory__
#define __RNAMake__pose_factory__

#include <stdio.h>

#include "motif/motif_factory.h"
#include "motif/pose.h"

class PoseFactory {
public:
    PoseFactory():
    mf_(MotifFactory())
    {}
    
    ~PoseFactory() {}
    
public:
    PoseOP
    pose_from_motif_tree(
        StructureOP const &,
        BasepairOPs const &,
        MotifOPs const &,
        std::map<Uuid, int, UuidCompare> const &);
    
private:
    void
    _add_motifs_to_pose(
        PoseOP &,
        MotifOPs const &);
    
    
private:
    MotifFactory mf_;
    
};


#endif /* defined(__RNAMake__pose_factory__) */
