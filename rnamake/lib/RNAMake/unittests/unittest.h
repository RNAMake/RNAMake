//
//  unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__unittest__
#define __RNAMake__unittest__

#include <stdexcept>
#include <stdio.h>
#include <iostream>
#include <map>

#include "util/settings.h"

String
unittest_resource_dir();

class UnittestException : public std::runtime_error {
public:
    UnittestException(
        String const & message) :
    std::runtime_error("Unittest Exception: " + message)
    {}
    
    int
    size() { return 0; }

};


class Unittest {
public:

public:
    
    virtual
    int run() { return  0; }
    
    virtual
    int run_all() { return 0; }
    
};

#endif /* defined(__RNAMake__unittest__) */
