//
//  motif_topology_unittests.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_topology.h"

#include "motif_topology_unittests.h"
#include "build/build_motif_graph.h"


namespace unittests {
namespace motif_structures {

void
MotifTopologyUnittest::test_to_tree() {
    auto builder = BuildMotifGraph();
    auto mg = builder.build(3);
    auto mt = graph_to_tree(mg);
}

int
MotifTopologyUnittest::run() {
    test_to_tree();
    return 1;
}

int
MotifTopologyUnittest::run_all() {
    return 1;
}

}
}
