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
namespace eternabot_unittests {

class EternabotStrategyUnittest : public Unittest {
public:
    
    EternabotStrategyUnittest() {}
    
    ~EternabotStrategyUnittest() {}
    
public:
    
    int
    test_creation();
    
    int
    test_pose_helix();
    
    int
    test_score();
    
    int
    test_score_compare();
    
public:
    
    int
    run();
    
    int
    run_all();
        
};


}
}

#endif /* defined(__RNAMake__eternabot_strategy_unittests__) */
