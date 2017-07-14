//
//  motif_state_search_solution.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_state.h"
#include "motif_state_search/motif_state_search_solution.h"


void
MotifStateSearchSolution::_get_path(
    MotifStateSearchNodeOP const & node) {
    
    auto current = node;
    while (current != nullptr) {
        path_.push_back(current);
        current = current->parent();
    }
    path_.pop_back();
    std::reverse(path_.begin(), path_.end());
}

MotifStateTreeOP
MotifStateSearchSolution::to_mst() {
    auto mst = std::make_shared<MotifStateTree>();
    mst->set_option_value("sterics", false);
    int i = 0, j = 0;
    for(auto const & n : path_) {
        auto cur_state = std::make_shared<MotifState>(n->ref_state()->name(),
                                                      n->ref_state()->end_names(),
                                                      n->ref_state()->end_ids(),
                                                      n->cur_state()->end_states(),
                                                      n->cur_state()->beads(),
                                                      n->ref_state()->score(),
                                                      n->ref_state()->size(),
                                                      n->ref_state()->block_end_add());
        
        // Nothing ever changes i...        
        if(i == 0) {
            mst->add_state(cur_state);
        }
        else {
            j = mst->add_state(cur_state, -1, n->parent_end_index());
            if(j == -1) {
                throw std::runtime_error("something went horribly wrong, cannot build solution");
            }
        }
        i++;
    }
    
    return mst;
}
