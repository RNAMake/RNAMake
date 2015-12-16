//
//  secondary_structure_tree_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_tree_unittest.h"
#include "build/build_secondary_structure.h"
#include "secondary_structure/secondary_structure_tree.h"



namespace unittests {
namespace sstruct_unittests  {

void
SecondaryStructureTreeUnittest::test_creation() {
    auto sst = sstruct::SecondaryStructureTree();
    
}
    
void
SecondaryStructureTreeUnittest::test_from_pose() {
    auto builder = BuildSecondaryStructure();
    auto p = builder.build_helix();
    auto sst = sstruct::tree_from_pose(p);
}

int
SecondaryStructureTreeUnittest::run() {
    test_creation();
    test_from_pose();
    return 0;
}



}
}