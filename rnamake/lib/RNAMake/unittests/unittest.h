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

class FailedTestException : public std::runtime_error {
public:
    FailedTestException(
        String const & message) :
    std::runtime_error("FailedTestException Exception: " + message)
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
    std::function<void()> f,
    String const & message) {
    
    try {
        f();
        throw FailedTestException("failed exception");
    }
    catch (T const & e) {}
    
    catch(FailedTestException const & e) {
        throw UnittestException(message);
    }
    
    catch(std::exception const & e) {
        throw UnittestException("unexcepted exception");
    }
    
}







#endif /* defined(__RNAMake__unittest__) */
