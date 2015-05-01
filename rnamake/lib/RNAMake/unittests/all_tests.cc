//
//  all_tests.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>


#include <vector>

#include "util_unittests/uuid_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_type_set_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/chain_unittest.h"
#include "structure_unittests/structure_unittest.h"

#include "resources_unittests/motif_library_unittest.h"
#include "resources_unittests/library_manager_unittest.h"

int run_util_unittests() {
    
    UuidUnittest uuid_test;
    uuid_test.run();
    
    return 1;

}


int run_structure_unittests() {
    
    std::vector<Unittest*> unittests(5);
    unittests[0] = new AtomUnittest();
    unittests[1] = new ResidueTypeSetUnittest();
    unittests[2] = new ResidueUnittest();
    unittests[3] = new ChainUnittest();
    unittests[4] = new StructureUnittest();
    
    for(auto & test : unittests) { test->run(); }
    
    return 1;
}

int run_resource_unittests() {
    
    std::vector<Unittest*> unittests(2);
    unittests[0] = new MotifLibraryUnittest();
    unittests[1] = new LibraryManagerUnittest();
    

}



int main(int argc, const char * argv[]) {
    
    run_util_unittests();
    run_structure_unittests();
    run_resource_unittests();
    return 0;
}