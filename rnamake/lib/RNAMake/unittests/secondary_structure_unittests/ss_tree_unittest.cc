//
//  ss_tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "ss_tree_unittest.h"
#include "secondary_structure/ss_tree.h"
int
SS_TreeUnittest::test_creation() {
    //SS_Tree ss_tree("((()))","GGGCCC");
    //std::cout << ss_tree.size() << std::endl;
    return 1;
}


int
SS_TreeUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    return 1;
}