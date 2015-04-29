//
//  motif_tree_state_path_refiner.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_path_refiner__
#define __RNAMake__motif_tree_state_path_refiner__

#include <stdio.h>
#include "motif_tree_state_search.h"


class MotifTreeStatePathRefiner {
public:
    MotifTreeStatePathRefiner():
    search_ (  MotifTreeStateSearch() ),
    options_ ( MotifTreeStateSearchOptions() ) {
        options_.numeric("max_size", 100);
        options_.numeric("max_node_level", 12);
        options_.numeric("max_n_solutions", 1);
    }
    
    ~MotifTreeStatePathRefiner() {}
    
public:
    MotifTreeStateSearchSolutionOPs
    find_path(
        BasepairStateOP const &,
        BasepairStateOP const &);
    
    void
    set_search_options();
    
public:
    
    inline
    void
    set_numeric_option(
        String const option,
        float const value) {
        options_.numeric(option, value);
    }
    
private:
    MotifTreeStateSearch search_;
    MotifTreeStateSearchOptions options_;
    
};


#endif /* defined(__RNAMake__motif_tree_state_path_refiner__) */
