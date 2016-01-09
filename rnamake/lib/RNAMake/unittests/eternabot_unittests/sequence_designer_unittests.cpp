//
//  sequence_designer_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "sequence_designer_unittests.h"

#include "build/build_secondary_structure.h"
#include "eternabot/sequence_designer.h"


namespace unittests {
namespace eternabot_unittests {


int
SequenceDesignerUnittest::test_creation() {
    auto designer = eternabot::SequenceDesigner();
    
    return 0;
}
    
int
SequenceDesignerUnittest::test_design() {
    auto designer = eternabot::SequenceDesigner();
    designer.set_option_value("steps", 1000);
    designer.setup();
    auto builder = BuildSecondaryStructure();
    auto p = builder.build_hairpin(20, 1);
    std::cout << p->sequence() << std::endl;
    std::cout << p->dot_bracket() << std::endl;
    auto results = designer.design(p);
    std::cout << results[0]->score << " " << results[1]->sequence << std::endl;
    return 0;
}
    
int
SequenceDesignerUnittest::run() {
    test_creation();
    test_design();
    return 0;
}

int
SequenceDesignerUnittest::run_all() {
    return 0;
}


    
}
}