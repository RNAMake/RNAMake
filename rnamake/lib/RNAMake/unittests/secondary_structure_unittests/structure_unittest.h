//
//  structure_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sec_structure_unittest__
#define __RNAMake__sec_structure_unittest__


#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace sstruct_unittests  {
    
class StructureUnittest : public Unittest {
public:
    StructureUnittest() {}
    
    ~StructureUnittest() {}
    
public:
    
    void
    test_creation();
    
    void
    test_find_residue();
    
    void
    test_copy();
    
    void
    test_to_str();
    
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};
    
}
}

#endif /* defined(__RNAMake__sec_structure_unittest__) */
