//
//  motif_assembly_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_assembly_unittest.h"
#include "motif_assembly/motif_assembly.h"

int
MotifAssemblyUnittest::test_add_motif() {
    MotifAssembly ma;
    
    return 1;
}


int
MotifAssemblyUnittest::run() {
    if (test_add_motif() == 0)    {  std::cout << "test_add_motif failed" << std::endl; }
    return 1;
}