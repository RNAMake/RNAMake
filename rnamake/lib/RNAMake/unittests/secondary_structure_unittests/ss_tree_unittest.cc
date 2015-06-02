//
//  ss_tree_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "ss_tree_unittest.h"
#include "secondary_structure/ss_tree.h"
#include "secondary_structure/ss_tree_node.h"


int
SS_TreeUnittest::test_creation() {
    //SS_Tree ss_tree("((()))","GGGCCC");
    //std::cout << ss_tree.size() << std::endl;
    
    //SS_Tree ss_tree2("(....)","GGGCCC");
    //std::cout << ss_tree2.size() << std::endl;

    String seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA";
    String ss  = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...)))))..";
    
    //SS_Tree ss_tree3(ss, seq);
    //std::cout << ss_tree3.size() << std::endl;

    //SS_Tree ss_tree4("..(.)..", "AACGGAA");
    
    SS_Tree ss_tree5("((...(.).(.)...))", "GGAAAGACAGACAAACC");
    //for(auto const & n : ss_tree5) {
    //    std::cout << n.index() << " " << n.data()->what() << " " << n.data()->seq() << std::endl;
    //}

    
    //std::cout << ss_tree.get_node(0)->data()->yb() << std::endl;
    
    return 1;
}


int
SS_TreeUnittest::run() {
    if (test_creation() == 0)  {  std::cout << "test_creation failed" << std::endl; }
    return 1;
}