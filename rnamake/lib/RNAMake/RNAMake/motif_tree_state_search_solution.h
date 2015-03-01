//
//  motif_tree_state_search_solution.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_search_solution__
#define __RNAMake__motif_tree_state_search_solution__

#include <stdio.h>
#include "motif_tree_state_search_node.h"
#include "motif_tree_state_tree.h"

class MotifTreeStateSearchSolution {
public:
    MotifTreeStateSearchSolution(
        MotifTreeStateSearchNodeOP const &,
        float);
    
public:
    
    MotifTreeStateTree
    to_mtst();
    
    
public:
    
    inline
    MotifTreeStateSearchNodeOPs const &
    path() const { return path_; }
    
    
private:
    MotifTreeStateSearchNodeOPs path_;
    float score_;
    
};

typedef std::shared_ptr<MotifTreeStateSearchSolution> MotifTreeStateSearchSolutionOP;
typedef std::vector<MotifTreeStateSearchSolutionOP> MotifTreeStateSearchSolutionOPs;

#endif /* defined(__RNAMake__motif_tree_state_search_solution__) */
