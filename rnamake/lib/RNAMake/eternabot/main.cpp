//
//  main.cpp
//  eternabot
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "vienna.h"
#include "sequence_designer.h"

int
test_design_sequence() {
    SequenceDesigner designer;
    designer.design("NNNNNUUCGNNNNN", "(((((....)))))");
    return 1;
}


int main(int argc, const char * argv[]) {
    // insert code here...
    
    test_design_sequence();
//    std::cout << "out" << std::endl;
    
    /*Vienna v;
    FoldResult fr = v.vfold("GGGAAACCC");
    std::cout << fr.structure << std::endl;
    FoldResult fr2 = v.vfold("GGGGGGGGGGGGGGGCCCCCCCCCCCCC");
    std::cout << fr2.free_energy << std::endl;*/


    return 0;
}
