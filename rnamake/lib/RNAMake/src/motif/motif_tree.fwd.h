//
//  motif_tree.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_motif_tree_fwd_h
#define RNAMake_motif_tree_fwd_h

#include "data_structure/graph/graph.h"
#include "data_structure/graph/graph_node.h"

typedef GraphNodeOP<MotifOP> MotifTreeNodeOP;
typedef GraphConnectionOP<MotifOP> MotifTreeConnectionOP;

class MotifTreeException : public std::runtime_error {
public:
    MotifTreeException(
        String const & message) :
    std::runtime_error(message)
    {}
    
};

#endif
