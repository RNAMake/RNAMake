//
//  tool_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/18/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "tool_unittests.hpp"


namespace unittests {
    
class TestException : public std::runtime_error {
public:
    TestException(
        String const & message):
    std::runtime_error("Test Exception: " + message)
    {}
};
    
    
void
ToolUnittest::test_failUnlessThrows() {
    std::function<void()> f = []() {  throw TestException("test"); };
    
    //throws correct error
    failUnlessThrows<TestException>(f, "test");
    
    //does not throw and error
    std::function<void()> f2 = []() { };
    
    try {
        failUnlessThrows<TestException>(f2, "test");
    }
    catch(UnittestException const & e) {}
    catch(...) {
        throw UnittestException("failUnlessThrows did not throw correct exception");
    }
    
    //throws wrong error
    std::function<void()> f3 = []() { throw UnittestException("test");  };
    
    try {
        failUnlessThrows<TestException>(f3, "test");
    }
    catch(UnittestException const & e) {}
    catch(...) {
        throw UnittestException("failUnlessThrows did not throw correct exception");
    }

}

int
ToolUnittest::run() {
    test_failUnlessThrows();
    return 0;
}
 
int
ToolUnittest::run_all() {
    return 0;

}

    
}