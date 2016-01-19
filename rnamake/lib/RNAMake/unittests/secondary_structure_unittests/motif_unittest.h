//
//  motif_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sec_motif_unittest__
#define __RNAMake__sec_motif_unittest__

#include <stdio.h>


#include "unittest.h"

namespace unittests {
namespace sstruct_unittests  {
    
class MotifUnittest : public Unittest {
public:
    MotifUnittest() {}
    
    ~MotifUnittest() {}
    
public:
    
    void
    test_creation();
    
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
    
} // sstruct_unittests
} // unittests


#endif /* defined(__RNAMake__motif_unittest__) */
