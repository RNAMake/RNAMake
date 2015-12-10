//
//  cl_option_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__cl_option_unittest__
#define __RNAMake__cl_option_unittest__

#include <stdio.h>

#include "unittest.h"

class CL_OptionUnittest : public Unittest {
public:
    CL_OptionUnittest() {}
    
    ~CL_OptionUnittest() {}
    
    int
    size() { return 3; }
    
public:
    
    int
    test_add_option();
    
    int
    test_parse_1();
    
    int
    test_parse_2();
    

public:
    
    int
    run();
    
    //void
    //run_all();
    
};

#endif /* defined(__RNAMake__cl_option_unittest__) */
