//
//  pose.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__pose__
#define __RNAMake__pose__

#include <map>
#include <stdio.h>
#include <memory>
#include "motif.h"

class Pose : public Motif {
public:
    Pose():
    designable_ (std::map<Uuid, int>()) {}
    Pose(MotifOP const &);
    Pose(StructureOP const &, BasepairOPs const &);

public:
    
    void
    set_bp_designable(BasepairOP const &);

public:
    
    inline
    void
    designable(std::map<Uuid, int> const & ndesignable) { designable_ = ndesignable; }
    
private:
    std::map<Uuid, int> designable_;

    
};

typedef std::shared_ptr<Pose> PoseOP;

#endif /* defined(__RNAMake__pose__) */

