//
//  basepair_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__basepair_unittest__
#define __RNAMake__basepair_unittest__

#include <stdio.h>

//RNAMake Headers
#include "unittest.h"
#include "structure/structure.h"
#include "structure/basepair.h"


class BasepairUnittest : public Unittest {
public:
    
    BasepairUnittest();
    
    ~BasepairUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_str_to_basepair_state();
    
    int
    test_basepair_state_to_str();
    
    int
    test_get_transforming_r_and_t_test();
    
    int
    test_move();
    
public:
    
    int
    run();
    
private:
    
    Structure s_;
 
    
};

#endif /* defined(__RNAMake__basepair_unittest__) */
