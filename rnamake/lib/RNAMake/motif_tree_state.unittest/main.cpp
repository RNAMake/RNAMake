//
//  main.cpp
//  motif_tree_state.unittest
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "motif_tree_state.h"
#include "motif_tree_state_library.h"
#include "motif_tree_state_node.h"
#include "motif_tree_state_node_aligner.h"

int
test_creation() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    return 1;
}

int
test_creation_node() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNode mtn ( mts_lib.motif_tree_states()[0], 0, NULL, 0);
    BasepairStateOPs avail_ends = mtn.available_ends();
    return 1;
}

int
test_add_child() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNodeOP mtn1 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, NULL, 0));
    MotifTreeStateNodeOP mtn2 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, mtn1, 0));
    BasepairStateOPs ends = mtn1->available_ends();
    mtn1->add_child(mtn2, ends[0]);
    ends = mtn1->available_ends();
    int parent_end_index = mtn2->parent_end_index();
    BasepairStateOP parent_end = mtn2->parent_end();
    return 1;
}

int
test_replace_mts() {
    MotifTreeStateLibrary mts_lib ( TWOWAY );
    MotifTreeStateNodeOP mtn1 (new MotifTreeStateNode( mts_lib.motif_tree_states()[0], 0, NULL, 0));
    mtn1->replace_mts(mts_lib.motif_tree_states()[1]);
    return 1;
}


int main(int argc, const char * argv[]) {
    if (test_creation() == 0)      { std::cout << "test_creation failed" << std::endl;  }
    if (test_creation_node() == 0) { std::cout << "test_creation_node failed" << std::endl;  }
    if (test_add_child() == 0)     { std::cout << "test_add_child failed" << std::endl;  }
    if (test_replace_mts() == 0)   { std::cout << "test_add_child failed" << std::endl;  }

    return 0;
}










