//
//  structure_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure_factory__
#define __RNAMake__structure_factory__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "structure/structure.h"
#include "structure/pdb_parser.h"

class StructureFactory {
    
public:
    
    StructureOP
    get_structure(
        String const & );
    
    
    ChainOPs
    build_chains(
        ResidueOPs & );
      
};

#endif /* defined(__RNAMake__structure_factory__) */
