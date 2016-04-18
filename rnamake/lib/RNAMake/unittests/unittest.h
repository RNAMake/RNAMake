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
#include <functional>


#include "util/settings.h"

String
unittest_resource_dir();

class UnittestException : public std::runtime_error {
public:
    UnittestException(
        String const & message) :
    std::runtime_error("Unittest Exception: " + message)
    {}
};


class Unittest {
public:

public:
    
    virtual
    int
    size() { return 0; }
    
    virtual
    int run() { return  0; }
    
    virtual
    int run_all() { return 0; }
    
    
};


inline
void
failUnless(
    bool statement,
    String const & message) {
    
    if(!statement) { throw UnittestException(message); }
    
}

template <typename T>
void
failUnlessThrows(
    std::function<void()> f) {
    
    try {
        
    }
    catch (T const & e) {
    }
    
    catch(...) {
        
    }
    
}







#endif /* defined(__RNAMake__unittest__) */
