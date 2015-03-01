//
//  main.cpp
//  secondary_structure_tree.unittest
//
//  Created by Joseph Yesselman on 2/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "secondary_structure_tree.h"

int
test_sstree() {
    SecondaryStructureTree sstree("((().))", "GGGCACC");
    SSandSeqOP ss_and_seq = sstree.get_ss_and_seq();
    String seq = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA";
    String ss =
        "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...)))))..";
    SecondaryStructureTree sstree2(ss, seq);
    ss_and_seq = sstree2.get_ss_and_seq();

    if(seq.compare(ss_and_seq->seq) != 0) {
        std::cout << "fail" << std::endl;
    }
    
    return 1;
}

int main(int argc, const char * argv[]) {
    if (test_sstree() == 0)          { std::cout << "test_sstree failed" << std::endl;  }
    return 0;
}










