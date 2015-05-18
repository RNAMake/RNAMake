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
#include "base_unittests/option_unittest.h"

#include "util_unittests/uuid_unittest.h"
#include "util_unittests/x3dna_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_type_set_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/chain_unittest.h"
#include "structure_unittests/structure_unittest.h"
#include "structure_unittests/basepair_unittest.h"
#include "structure_unittests/pdb_parser_unittest.h"

#include "motif_unittests/motif_unittest.h"
#include "motif_unittests/motif_tree_unittest.h"
#include "motif_unittests/motif_tree_merger_unittest.h"
#include "motif_unittests/motif_scorer_unittest.h"

#include "motif_tree_state_unittests/motif_tree_state_library_unittest.h"

#include "resources_unittests/motif_library_unittest.h"
#include "resources_unittests/library_manager_unittest.h"


int main(int argc, const char * argv[]) {
    std::vector<Unittest*> units;
    units.push_back(new OptionUnittest());
    units.push_back(new UuidUnittest());
    units.push_back(new X3dnaUnittest());
    units.push_back(new AtomUnittest());
    units.push_back(new ResidueTypeSetUnittest());
    units.push_back(new ResidueUnittest());
    units.push_back(new ChainUnittest());
    units.push_back(new StructureUnittest());
    units.push_back(new BasepairUnittest());
    units.push_back(new PDBParserUnittest());
    units.push_back(new MotifUnittest());
    units.push_back(new MotifTreeUnittest());
    units.push_back(new MotifTreeMergerUnittest());
    units.push_back(new MotifScorerUnittest());
    units.push_back(new MotifTreeStateLibraryUnittest());
    units.push_back(new MotifLibraryUnittest());
    units.push_back(new LibraryManagerUnittest());
    
    for(auto const & test : units) {
        if(test == nullptr) { continue; }
        test->run_all();
    }

    
    
    return 0;
}












