//
//  tool_unittests.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/18/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef tool_unittests_hpp
#define tool_unittests_hpp

#include <stdio.h>

#include "unittest.h"

namespace unittests {
    
class ToolUnittest : public Unittest {
public:
    
    ToolUnittest() {}
    
    ~ToolUnittest() {}
    
public:
    
    void
    test_failUnlessThrows();
    
public:
    
    int
    run();
    
    int
    run_all();
    
private:
    
};
    
}

#endif /* tool_unittests_hpp */
