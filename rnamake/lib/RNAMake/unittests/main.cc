//
//  main.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "base_unittests/option_unittest.h"
#include "base_unittests/cl_option_unittest.h"

#include "data_structure_unittests/tree_unittest.h"
#include "data_structure_unittests/graph_unittest.h"

#include "util_unittests/sqlite3_connection_unittest.h"
#include "util_unittests/x3dna_unittest.h"
#include "util_unittests/uuid_unittest.h"

#include "secondary_structure_unittests/secondary_structure_factory_unittest.h"
#include "secondary_structure_unittests/secondary_structure_unittest.h"
#include "secondary_structure_unittests/ss_tree_unittest.h"

#include "structure_unittests/structure_unittest.h"

#include "motif_unittests/motif_unittest.h"
#include "motif_unittests/motif_factory_unittest.h"
#include "motif_unittests/motif_tree_unittest.h"
#include "motif_unittests/motif_tree_merger_unittest.h"
#include "motif_unittests/pose_factory_unittest.h"
#include "motif_unittests/motif_state_unittest.h"

#include "resources_unittests/motif_sqlite_connection_unittest.h"
#include "resources_unittests/motif_sqlite_library_unittest.h"
#include "resources_unittests/motif_state_sqlite_library_unittest.h"
#include "resources_unittests/motif_state_ensemble_sqlite_library_unittest.h"
#include "resources_unittests/resource_manager_unittest.h"
#include "resources_unittests/added_motif_library_unittest.h"

#include "motif_data_structures_unittests/motif_state_tree_unittest.h"
#include "motif_data_structures_unittests/motif_state_ensemble_unittest.h"

/*#include "util_unittests/uuid_unittest.h"


#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/resource_manager_unittest.h"
#include "structure_unittests/pdb_parser_unittest.h"
#include "structure_unittests/basepair_unittest.h"

#include "motif_unittests/motif_tree_merger_unittest.h"
#include "motif_unittests/motif_scorer_unittest.h"
#include "motif_unittests/motif_tree_topology_unittest.h"

#include "motif_tree_state_unittests/motif_tree_state_library_unittest.h"

#include "resources_unittests/motif_library_unittest.h"
#include "resources_unittests/library_manager_unittest.h"

#include "motif_assembly_unittests/motif_assembly_unittest.h"*/


int main(int argc, const char * argv[]) {
    
    MotifStateEnsembleUnittest test;
    test.run();
}