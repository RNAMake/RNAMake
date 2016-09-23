//
//  instance_unittest.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/13/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef instance_unittest_hpp
#define instance_unittest_hpp

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace instances {

class InstanceUnittests : public Unittest {
public:
    InstanceUnittests() {}
    
    ~InstanceUnittests() {}
    
public:
    
    void
    test_residue();
    
public:
    
    int
    run();
    
    int
    run_all();
    
};
    
    
}
}

#endif /* instance_unittest_hpp */
