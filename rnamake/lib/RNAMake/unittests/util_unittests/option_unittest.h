//
//  option_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__option_unittest__
#define __RNAMake__option_unittest__

#include <stdio.h>

#include "unittest.h"

class OptionUnittest : public Unittest {
public:
    OptionUnittest() {}
    
    ~OptionUnittest() {}
    
public:
    
    int
    test_creation();

    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};


#endif /* defined(__RNAMake__option_unittest__) */
