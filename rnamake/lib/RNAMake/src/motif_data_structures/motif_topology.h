//
//  motif_toplogy.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_topology__
#define __RNAMake__motif_topology__

#include <stdio.h>

#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_graph.h"

MotifTreeOP
graph_to_tree(
    MotifGraphOP const & mg,
    GraphNodeOP<MotifOP> start = nullptr,
    String last_end = "");


#endif /* defined(__RNAMake__motif_toplogy__) */
