//
//  motif_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_factory__
#define __RNAMake__motif_factory__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "util/settings.h"
#include "structure/structure_factory.h"
#include "motif/motif_to_secondary_structure.h"
#include "motif/motif_scorer.h"
#include "motif/motif.h"


class MotifFactory {
public:
    MotifFactory():
    sf_(StructureFactory()),
    parser_(MotiftoSecondaryStructure()) {
        auto path = motif_dirs() + "ref.motif";
        ref_motif_ = file_to_motif(path);
        path = motif_dirs() + "base.motif";
        base_motif_ = file_to_motif(path);
        added_helix_ = std::make_shared<Motif>(base_motif_->copy());
        clash_radius_ = 2.9;
    }
    
    ~MotifFactory() {}
    
public:
    
    MotifOP
    motif_from_file(
        String const & path);
    
    void
    standardize_motif(
        MotifOP &);
    
    MotifOP
    can_align_motif_to_end(
        MotifOP const &,
        int);
    
    MotifOP
    align_motif_to_common_frame(
        MotifOP const &,
        int);
    
private:
    
    BasepairOPs
    _setup_basepairs(
        String const &,
        StructureOP const &);
    
    BasepairOPs
    _setup_basepair_ends(
        StructureOP const &,
        BasepairOPs const &);
    
    void
    _setup_secondary_structure(
        MotifOP &);
    
    void
    _align_chains(
        MotifOP &);
    
    void
    _align_ends(
        MotifOP &);
    
    int
    _steric_clash(
        MotifOP const &,
        MotifOP const &);
    
private:
    StructureFactory sf_;
    MotiftoSecondaryStructure parser_;
    MotifScorer scorer_;
    MotifOP ref_motif_, base_motif_, added_helix_;
    float clash_radius_;
    
    
};

#endif /* defined(__RNAMake__motif_factory__) */
