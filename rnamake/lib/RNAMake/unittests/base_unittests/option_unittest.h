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

namespace unittests {
    
class OptionUnittest : public Unittest {
public:
    OptionUnittest() {}
    
    ~OptionUnittest() {}
    
    int
    size() { return 3; }
    
public:
    
    int
    test_creation();
    
    int
    test_add_option();
    
    int
    test_option();
    
public:
    
    int
    run();
    
    int
    run_all();
    
};
    
}


#endif /* defined(__RNAMake__option_unittest__) */
