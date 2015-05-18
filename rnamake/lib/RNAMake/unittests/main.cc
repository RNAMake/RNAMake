//
//  main.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util_unittests/uuid_unittest.h"
#include "util_unittests/x3dna_unittest.h"
#include "util_unittests/option_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/resource_manager_unittest.h"
#include "structure_unittests/pdb_parser_unittest.h"
#include "structure_unittests/structure_unittest.h"
#include "structure_unittests/basepair_unittest.h"

#include "motif_unittests/motif_unittest.h"
#include "motif_unittests/motif_tree_unittest.h"
#include "motif_unittests/motif_tree_merger_unittest.h"
#include "motif_unittests/motif_scorer_unittest.h"

#include "motif_tree_state_unittests/motif_tree_state_library_unittest.h"

#include "resources_unittests/motif_library_unittest.h"
#include "resources_unittests/library_manager_unittest.h"


int main(int argc, const char * argv[]) {

    OptionUnittest test;
    test.run();
    

}