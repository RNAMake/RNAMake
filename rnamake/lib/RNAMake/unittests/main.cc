//
//  main.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <vector>

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_type_set_unittest.h"

int run_structure_unittests() {
    
    std::vector<Unittest*> unittests(2);
    unittests[0] = new AtomUnittest();
    unittests[1] = new ResidueTypeSetUnittest();
    
    for(auto & test : unittests) { test->run(); }
    
    return 1;
}



int main(int argc, const char * argv[]) {

    run_structure_unittests();


}