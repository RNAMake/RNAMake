//
//  path_follower_unittests.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__path_follower_unittests__
#define __RNAMake__path_follower_unittests__

#include <stdio.h>

#include "unittest.h"

namespace unittests {
namespace motif_state_search {

class PathFollowerUnittest : public Unittest {
public:
    PathFollowerUnittest() {}
    
    ~PathFollowerUnittest() {}
    
public:
    
    int
    test_creation();

    int
    test_build();
    
public:
    
    int
    run();
    
    //void
    //run_all();
    
};
        
}
}

#endif /* defined(__RNAMake__path_follower_unittests__) */
