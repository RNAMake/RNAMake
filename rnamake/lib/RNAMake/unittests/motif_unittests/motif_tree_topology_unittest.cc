//
//  motif_tree_topology_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_topology_unittest.h"
#include "motif/motif_tree_topology.h"

int
MotifTreeTopologyUnittest::test_creation() {
    SS_Tree ss_t("(((+)))", "GAG+CUC");
    MotifTreeTopology mtt(ss_t);
    
    return 1;
}


int
MotifTreeTopologyUnittest::run() {
    if (test_creation() == 0)            { std::cout << "test_creation failed" << std::endl;  }
    return 0;
}
