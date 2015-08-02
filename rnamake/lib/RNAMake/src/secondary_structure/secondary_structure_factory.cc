//
//  secondary_structure_factory.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_factory.h"

namespace sstruct {

void
SecondaryStructureFactory::_get_basepairs(
    SS_Tree const & sstree,
    SecondaryStructureOP & ss) {
    
    for(auto const & n : sstree) {
        if(n->data()->type() != SS_NodeData::SS_Type::SS_BP) { continue; }
        
        
        
    }
    
}


} //sstruct