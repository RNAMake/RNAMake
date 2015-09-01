//
//  motif_factory_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/9/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_factory_unittest__
#define __RNAMake__motif_factory_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"

class MotifFactoryUnittest : public Unittest {
public:
    
    MotifFactoryUnittest() {}
    
    ~MotifFactoryUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_load();
    
    int
    test_standardize_motif();
    
    int
    test_can_align_motif_to_end();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__motif_factory_unittest__) */
