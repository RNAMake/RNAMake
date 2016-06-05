//
//  util.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__util__
#define __RNAMake__util__

#include <stdio.h>


//RNAMake Headers
#include "secondary_structure/rna_structure.h"
#include "secondary_structure/pose.h"

namespace sstruct {


String
assign_end_id(
    RNAStructureOP const &,
    BasepairOP const &);


void
fill_basepairs_in_ss(PoseOP &);
    
    
}
#endif /* defined(__RNAMake__util__) */
