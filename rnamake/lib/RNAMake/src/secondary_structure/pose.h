//
//  pose.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake_sec_pose__
#define __RNAMake_sec_pose__

#include <stdio.h>

// RNAMake Headers
#include "secondary_structure/rna_structure.h"
#include "secondary_structure/motif.h"

namespace sstruct {

class Pose : public RNAStructure {
public:
    Pose():
    RNAStructure()
    {}
    
    Pose(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    RNAStructure(structure, basepairs, ends)
    {}
    
    Pose(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends,
        MotifOPs const & motifs):
    RNAStructure(structure, basepairs, ends),
    motifs_(motifs)
    {}
    
    Pose(
        RNAStructureOP const & rna_struc,
        MotifOPs const & motifs):
    RNAStructure(*rna_struc),
    motifs_(motifs)
    {}
    
public:
    MotifOPs const &
    motifs() { return motifs_; }

private:
    MotifOPs motifs_;
    
};

typedef std::shared_ptr<Pose> PoseOP;
    
}
#endif /* defined(__RNAMake__pose__) */
