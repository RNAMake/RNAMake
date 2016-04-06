//
//  application_unittest.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/1/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef application_unittest_hpp
#define application_unittest_hpp

#include <stdio.h>


#include "unittest.h"

namespace unittests {
    
class ApplicationUnittest : public Unittest {
public:
    ApplicationUnittest() {}
    
    ~ApplicationUnittest() {}
    
    int
    size() { return 4; }
    
public:
    
    int
    test_creation();
    
    int
    test_setup_cl_options();
    
public:
    
    int
    run();
    
    int
    run_all();
    
};
    
}


#endif /* application_unittest_hpp */
