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
#include "structure/structure.h"
#include "structure/pdb_parser.h"


class MotifFactory {
public:
    MotifFactory() {}
    
    ~MotifFactory() {}
    
public:
    
    StructureOP
    get_structure(
        String const & );
    
private:
    
    ChainOPs
    _build_chains(
        ResidueOPs & );
    
    
};

#endif /* defined(__RNAMake__motif_factory__) */
