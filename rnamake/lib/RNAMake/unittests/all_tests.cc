//
//  all_tests.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <stdio.h>
#include <vector>

//RNAMake Headers
#include "util_unittests/uuid_unittest.h"
#include "util_unittests/x3dna_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_type_set_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/chain_unittest.h"
#include "structure_unittests/structure_unittest.h"
#include "structure_unittests/basepair_unittest.h"

#include "resources_unittests/motif_library_unittest.h"
#include "resources_unittests/library_manager_unittest.h"

int run_resource_unittests() {
    
    std::vector<Unittest*> unittests(2);
    unittests[0] = new MotifLibraryUnittest();
    unittests[1] = new LibraryManagerUnittest();
    return 1;
    

}

int main(int argc, const char * argv[]) {
    std::vector<Unittest*> units;
    units.push_back(new UuidUnittest());
    units.push_back(new X3dnaUnittest());
    units.push_back(new AtomUnittest());
    units.push_back(new ResidueTypeSetUnittest());
    units.push_back(new ResidueUnittest());
    units.push_back(new ChainUnittest());
    units.push_back(new StructureUnittest());
    units.push_back(new BasepairUnittest());

    for(auto const & test : units) {
        if(test == nullptr) { continue; }
        test->run_all();
    }

    
    
    return 0;
}












