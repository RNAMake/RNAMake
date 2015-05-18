//
//  motif_scorer_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_scorer_unittest__
#define __RNAMake__motif_scorer_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifScorerUnittest : public Unittest {
public:
    
    MotifScorerUnittest() {}
    
    ~MotifScorerUnittest() {}
    
public:
    
    int
    test_score();
    
public:
    
    int
    run();
    
    void
    run_all();
    
};

#endif /* defined(__RNAMake__motif_scorer_unittest__) */
