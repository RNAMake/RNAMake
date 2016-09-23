//
//  secondary_structure_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__secondary_structure_tree_unittest__
#define __RNAMake__secondary_structure_tree_unittest__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace sstruct_unittests  {

class SecondaryStructureTreeUnittest : public Unittest {
public:
    SecondaryStructureTreeUnittest() {}
    
    ~SecondaryStructureTreeUnittest() {}
    
public:
    
    void
    test_creation();
    
    void
    test_from_pose();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};
    
}
}


#endif /* defined(__RNAMake__secondary_structure_tree_unittest__) */
