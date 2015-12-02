//
//  motif_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_unittest.h"
#include "secondary_structure/secondary_structure_factory.h"

namespace unittests {
namespace sstruct_unittests  {

void
MotifUnittest::test_creation() {
    sstruct::SecondaryStructureFactory ssf;
    auto m = ssf.motif("GGGGG+CCCCC", "(((((+)))))");
}
    

int
MotifUnittest::run() {
    test_creation();
    return 0;
}
    
} // sstruct
} // sstruct_unittests