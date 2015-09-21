//
//  motif_state_search_solution.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_solution__
#define __RNAMake__motif_state_search_solution__

#include <stdio.h>

//RNAMake Headers
#include "motif_data_structures/motif_state_tree.h"
#include "motif_state_search/motif_state_search_node.h"

class MotifStateSearchSolution {
public:
    MotifStateSearchSolution(
        MotifStateSearchNodeOP const & node,
        float score):
    score_(score) {
        _get_path(node);
    }
    
    ~MotifStateSearchSolution() {}
    
public:
    
    MotifStateTreeOP
    to_mst();
    
private:
    
    void
    _get_path(
        MotifStateSearchNodeOP const &);

private:
    MotifStateSearchNodeOPs path_;
    float score_;

};

typedef std::shared_ptr<MotifStateSearchSolution> MotifStateSearchSolutionOP;
typedef std::vector<MotifStateSearchSolutionOP>   MotifStateSearchSolutionOPs;


#endif /* defined(__RNAMake__motif_state_search_solution__) */
