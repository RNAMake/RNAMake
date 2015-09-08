//
//  motif_state_ensemble_tree_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_tree_unittest__
#define __RNAMake__motif_state_ensemble_tree_unittest__

#include <stdio.h>


//RNAMake Headers
#include "unittest.h"

class MotifStateEnsembleTreeUnittest : public Unittest {
public:
    
    MotifStateEnsembleTreeUnittest() {}
    
    ~MotifStateEnsembleTreeUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_add_ensemble();
    
    int
    test_from_mt();
    
    int
    test_to_mst();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__motif_state_ensemble_tree_unittest__) */
