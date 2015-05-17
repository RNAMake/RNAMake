//
//  unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__unittest__
#define __RNAMake__unittest__

#include <stdio.h>
#include <iostream>
#include <map>

#include "util/settings.h"

String
unittest_resource_dir();

class Unittest {
public:

public:
    
    virtual
    int run() { return  0; }
    
    virtual
    void run_all() { }
    
};

#endif /* defined(__RNAMake__unittest__) */
