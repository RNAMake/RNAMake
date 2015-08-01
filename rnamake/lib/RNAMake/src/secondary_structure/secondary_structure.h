//
//  secondary_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure__
#define __RNAMake__secondary_structure__

#include <stdio.h>

#include "base/types.h"
#include "util/uuid.h"

#include "secondary_structure/chain.h"
#include "secondary_structure/motif.h"

namespace sstruct {

    
class SecondaryStructure : public Motif {
public:
    
    SecondaryStructure() {}
    
    SecondaryStructure(
        String const &,
        String const &);
    
    SecondaryStructure(
        ChainOPs const &);
    
    ~SecondaryStructure() {}
    
public:
    
    
};
    
    
}



#endif /* defined(__RNAMake__secondary_structure__) */
