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
test_fold() {
    
    Vienna v;
    FoldResult fr = v.vfold("GGGAAACCC");
    //std::cout << fr.structure << std::endl;
    FoldResult fr2 = v.vfold("GGGGGGGGGGGGGGGCCCCCCCCCCCCC");
    //std::cout << fr2.free_energy << std::endl;
    return 1;
}

int
test_cofold() {
    Vienna v;
    FoldResult fr = v.vcofold("GGGAAACCCCCC&AAAAAGGGAAA");
    //std::cout << fr.structure << std::endl;
    return 1;
}

int
test_design_sequence() {
    SequenceDesigner designer;
    designer.design("NNNNNUUCGNNNNN", "(((((....)))))");
    return 1;
}

int
test_design_sequence_2() {
    SequenceDesigner designer;
    designer.design("NNNNNNNNNNNNNNNNNNNNNNNNN&NNNNNNNNNNNNNNNNNNNNNNNNN",
                    "(((((((((((((((((((((((((&)))))))))))))))))))))))))");
    return 1;
}


int main(int argc, const char * argv[]) {
    
    if (test_fold() == 0)             { std::cout << "test_fold failed" << std::endl;  }
    if (test_cofold() == 0)           { std::cout << "test_cofold failed" << std::endl;  }
    //if (test_design_sequence() == 0)  { std::cout << "test_design_sequence failed" << std::endl;  }
    if (test_design_sequence_2() == 0)  { std::cout << "test_design_sequence failed" << std::endl;  }

    
    return 0;
}
