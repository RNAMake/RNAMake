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

//RNAMake Headers
#include "motif/motif.h"

class Pose : public Motif {
public:
    Pose():
    designable_ (std::map<Uuid, int, UuidCompare>()) {}
    
    Pose(MotifOP const &);
    
    Pose(StructureOP const &, BasepairOPs const &);
    
    ~Pose() { }

public:
    
    void
    set_bp_designable(BasepairOP const &);
    
    String
    designable_sequence();

public:
    
    inline
    void
    designable(std::map<Uuid, int, UuidCompare> const & ndesignable) { designable_ = ndesignable; }
    
private:
    std::map<Uuid, int, UuidCompare> designable_;

    
};

typedef std::shared_ptr<Pose> PoseOP;

#endif /* defined(__RNAMake__pose__) */

