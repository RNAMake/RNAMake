//
//  path_follower_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//


#include "motif_state_search/path_follower.h"
#include "path_follower_unittests.h"


namespace unittests {
namespace motif_state_search {

int
PathFollowerUnittest::test_creation() {
    auto pf = PathFollower();
    return 0;
}
    
int
PathFollowerUnittest::test_build() {
    auto pf = PathFollower();
    
    return 0;
}
    
int
PathFollowerUnittest::run() {
    test_creation();
    test_build();
    return 0;
}
    

}
}