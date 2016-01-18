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

#include "secondary_structure_unittests/structure_unittest.h"
#include "secondary_structure_unittests/motif_unittest.h"
#include "secondary_structure_unittests/secondary_structure_parser_unittest.h"
#include "secondary_structure_unittests/secondary_structure_factory_unittest.h"
#include "secondary_structure_unittests/secondary_structure_tree_unittest.h"

#include "structure_unittests/atom_unittest.h"
#include "structure_unittests/residue_unittest.h"
#include "structure_unittests/chain_unittest.h"
#include "structure_unittests/structure_unittest.h"

#include "motif_unittests/motif_unittest.h"
#include "motif_unittests/motif_factory_unittest.h"

#include "motif_data_structures_unittests/motif_graph_unittest.h"
#include "motif_data_structures_unittests/motif_tree_unittest.h"
#include "motif_data_structures_unittests/motif_topology_unittests.h"

#include "motif_state_search_unittests/motif_state_search_unittest.h"
#include "motif_state_search_unittests/path_follower_unittests.h"

#include "resources_unittests/motif_sqlite_connection_unittest.h"
#include "resources_unittests/motif_sqlite_library_unittest.h"
#include "resources_unittests/resource_manager_unittest.h"

#include "eternabot_unittests/eternabot_strategy_unittests.h"
#include "eternabot_unittests/scorer_unittest.h"
#include "eternabot_unittests/sequence_designer_unittests.h"

#include "sequence_optimizer_unittests/sequence_optimizer_unittests.h"


int main(int argc, const char * argv[]) {
    
    unittests::motif_state_search::PathFollowerUnittest test;
    test.run();

    return 0;
}







