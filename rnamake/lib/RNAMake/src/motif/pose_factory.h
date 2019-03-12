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

#include "util/x3dna.h"
#include "motif/motif_factory.h"
#include "motif/pose.h"

class PoseFactory {
public:
    PoseFactory():
    mf_(MotifFactory()),
    clash_radius_(2.9)
    {}
    
    ~PoseFactory() {}
    
public:
    PoseOP
    pose_from_motif_tree(
        structure::StructureOP const &,
        structure::BasepairOPs const &,
        MotifOPs const &,
        std::map<util::Uuid, int, util::UuidCompare> const &);
    
    PoseOP
    pose_from_file(
        String const & path,
        int gu_are_helix=1,
        int signlet_bp_seperation=0);
    
private:
    void
    _setup_motifs_from_x3dna(
        PoseOP &,
        int,
        int);
    
    void
    _add_motifs_to_pose(
        PoseOP &,
        MotifOPs const &);
    
    void
    _add_secondary_structure_motifs(
        PoseOP &);
    
    void
    _standardize_prepose(
        MotifOP &);
    
    int
    _steric_clash(
        MotifOP const &,
        MotifOP const &);
    
    MotifOP
    _convert_x3dna_to_motif(
        util::X3Motif const &,
        PoseOP const &);
    
    
    
private:
    MotifFactory mf_;
    float clash_radius_;
    
};


#endif /* defined(__RNAMake__pose_factory__) */
