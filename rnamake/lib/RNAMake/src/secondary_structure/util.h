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

namespace sstruct {


String
assign_end_id(
    RNAStructureOP const &,
    BasepairOP const &);


}
#endif /* defined(__RNAMake__util__) */
