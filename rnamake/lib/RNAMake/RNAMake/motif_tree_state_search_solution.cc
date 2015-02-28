//
//  motif_tree_state_search_solution.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_state_search_solution.h"

MotifTreeStateSearchSolution::MotifTreeStateSearchSolution(
    MotifTreeStateSearchNodeOP const & node,
    float score):
    score_ (score )
{
    path_ = MotifTreeStateSearchNodeOPs();
    MotifTreeStateSearchNodeOP current = node;
    while ( current != NULL) {
        path_.push_back(current);
        current = current->parent();
    }
    
    std::reverse(path_.begin(), path_.end());
    
}

