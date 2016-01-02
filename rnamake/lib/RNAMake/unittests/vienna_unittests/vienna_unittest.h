//
//  vienna_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__vienna_unittest__
#define __RNAMake__vienna_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

namespace unittests {
namespace vienna {

class ViennaUnittest : public Unittest {
public:
    
    ViennaUnittest() {}
    
    ~ViennaUnittest() {}
    
    int
    size() { return 5; }
    
public:
    
    int
    test_creation();

    int
    test_folding();
    
    int
    test_folding_no_reset();
    
    int
    test_bp_probs();
    
    int
    test_memory_leak();
    
    
public:
    
    int
    run();
    
    int
    run_all();
    
};

}
}


#endif /* defined(__RNAMake__vienna_unittest__) */
