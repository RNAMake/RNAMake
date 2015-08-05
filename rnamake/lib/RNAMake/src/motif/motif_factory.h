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
#include "structure/structure_factory.h"
#include "motif/motif_to_secondary_structure.h"
#include "motif/motif_scorer.h"
#include "motif/motif.h"


class MotifFactory {
public:
    MotifFactory():
    sf_(StructureFactory()),
    parser_(MotiftoSecondaryStructure()){}
    
    ~MotifFactory() {}
    
public:
    
    MotifOP
    motif_from_file(
        String const & path);
    
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
    
    
private:
    StructureFactory sf_;
    MotiftoSecondaryStructure parser_;
    MotifScorer scorer_;
    
    
};

#endif /* defined(__RNAMake__motif_factory__) */
