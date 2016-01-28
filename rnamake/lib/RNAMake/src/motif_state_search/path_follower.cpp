//
//  path_follower.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/18/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_state_search/path_follower.h"

MotifTreeOP
PathFollower::next() {
    auto sol = search_.next();
    auto mt = sol->to_motif_tree();
    
    return mt;
    
}
