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
#include "motif/motif.h"


class MotifFactory {
public:
    MotifFactory():
    sf_(StructureFactory()) {}
    
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
    
    
private:
    StructureFactory sf_;
    
    
};

#endif /* defined(__RNAMake__motif_factory__) */
