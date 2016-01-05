//
//  scorer_unittest.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__eternabot_scorer_unittest__
#define __RNAMake__eternabot_scorer_unittest__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace eternabot_unittests {

class ScorerUnittest : public Unittest {
public:
    
    ScorerUnittest() {}
    
    ~ScorerUnittest() {}
    
public:
    
    int
    test_creation();
    
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


#endif /* defined(__RNAMake__scorer_unittest__) */
