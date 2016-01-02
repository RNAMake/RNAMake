//
//  eternabot_strategy_unittests.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__eternabot_strategy_unittests__
#define __RNAMake__eternabot_strategy_unittests__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace eternabot {

class EternabotStrategyUnittest : public Unittest {
public:
    
    EternabotStrategyUnittest() {}
    
    ~EternabotStrategyUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_score();
    
public:
    
    int
    run();
    
    int
    run_all();
        
};


}
}

#endif /* defined(__RNAMake__eternabot_strategy_unittests__) */
