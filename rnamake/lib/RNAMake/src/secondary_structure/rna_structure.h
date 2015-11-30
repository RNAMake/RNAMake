//
//  rna_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__rna_structure__
#define __RNAMake__rna_structure__

#include <stdio.h>

//RNAMake
#include "secondary_structure/basepair.h"
#include "secondary_structure/structure.h"

namespace sstruct {


class RNAStructure {
public:
    RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    structure_(structure),
    basepairs_(basepairs),
    ends_(end),
    name_(""),
    path_(""),
    score_(0),
    end_ids_(Strings())
    {}

private:
    StructureOP structure_;
    BasepairOPs basepairs_, ends_;
    String name_, path_;
    Strings end_ids_;
    float score_;
    
};

    
}

#endif /* defined(__RNAMake__rna_structure__) */
