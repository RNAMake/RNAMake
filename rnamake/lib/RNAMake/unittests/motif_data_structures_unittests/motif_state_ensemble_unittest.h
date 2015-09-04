//
//  motif_state_ensemble_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_unittest__
#define __RNAMake__motif_state_ensemble_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifStateEnsembleUnittest : public Unittest {
public:
    
    MotifStateEnsembleUnittest() {}
    
    ~MotifStateEnsembleUnittest() {}
    
public:
    
    int
    test_creation();

public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__motif_state_ensemble_unittest__) */
