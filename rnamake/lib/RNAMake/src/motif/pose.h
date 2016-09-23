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
#include "util/motif_type.h"
#include "motif/motif.h"

class Pose : public Motif {
public:
    Pose():
    designable_ (std::map<Uuid, int, UuidCompare>()),
    motifs_(std::map<MotifType, MotifOPs>())
    {}
    
    Pose(MotifOP const &);
    
    Pose(StructureOP const &, BasepairOPs const &);
    
    ~Pose() { }
    
public:
    
    MotifOPs const &
    motifs(
        MotifType const &);

public:
    
    void
    set_bp_designable(BasepairOP const &);
    
    String
    designable_sequence();

public:
    
    inline
    void
    designable(std::map<Uuid, int, UuidCompare> const & ndesignable) { designable_ = ndesignable; }
    
    inline
    void
    set_motifs(std::map<MotifType, MotifOPs> const & motifs) { motifs_ = motifs; }
    
    /*inline
    void
    set_ss_motifs(std::map<String, sstruct::MotifOPs> & ss_motifs) {
        secondary_structure_->set_motifs(ss_motifs);
    }*/
    
private:
    std::map<Uuid, int, UuidCompare> designable_;
    std::map<MotifType, MotifOPs> motifs_;

    
};

typedef std::shared_ptr<Pose> PoseOP;

#endif /* defined(__RNAMake__pose__) */

