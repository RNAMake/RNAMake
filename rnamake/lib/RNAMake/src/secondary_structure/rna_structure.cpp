//
//  rna_structure.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure/rna_structure.h"

namespace sstruct {

BasepairOPs
RNAStructure::get_basepair(
    Uuid const & bp_uuid) {
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->uuid() == bp_uuid) { bps.push_back(bp); }
        if(bp->res1()->uuid() == bp_uuid || bp->res2()->uuid() == bp_uuid ) { bps.push_back(bp); }
        
    }
    return bps;
}

BasepairOPs
RNAStructure::get_basepair(
    ResidueOP const & res1,
    ResidueOP const & res2) {
    BasepairOPs bps;
    for(auto & bp : basepairs_) {
        if(bp->res1()->uuid() == res1->uuid() && bp->res2()->uuid() == res2->uuid()) { bps.push_back(bp); }
        if(bp->res1()->uuid() == res2->uuid() && bp->res2()->uuid() == res1->uuid()) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
RNAStructure::get_basepair(
    Uuid const & uuid1,
    Uuid const & uuid2) {
    
    BasepairOPs bps;
    for( auto const & bp : basepairs_) {
        if(bp->res1()->uuid() == uuid1 && bp->res2()->uuid() == uuid2) { bps.push_back(bp); }
        if(bp->res1()->uuid() == uuid2 && bp->res2()->uuid() == uuid1) { bps.push_back(bp); }
    }
    return bps;
}

BasepairOPs
RNAStructure::get_basepair(
    String const & name) {
    for(auto const & bp : basepairs_) {
        if(name.compare(bp->name()) == 0) {
            return BasepairOPs{ bp };
        }
    }
    
    throw "could not find basepair with name " + name;
    
}
    
}
