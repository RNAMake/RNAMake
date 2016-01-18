//
//  sequence_optimizer_unittests.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_optimizer_unittests__
#define __RNAMake__sequence_optimizer_unittests__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace sequence_optimizer {

class SequenceOptimizerUnittest : public Unittest {
public:
    SequenceOptimizerUnittest() {}
    
    ~SequenceOptimizerUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_optimize();

public:
    
    int
    run();
    
    int
    run_all();
    
};
    
}
}
#endif /* defined(__RNAMake__sequence_optimizer_unittests__) */
