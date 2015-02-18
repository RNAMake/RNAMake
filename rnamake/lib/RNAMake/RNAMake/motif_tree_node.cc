//
//  motif_tree_node.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree_node.h"
#include "motif.h"
#include "types.h"
#include "settings.h"
#include "residue_type_set.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MotifTreeNode
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MotifTreeNode::MotifTreeNode(
    MotifOP const & m,
    int const level,
    int const index,
    int const flip):
    motif_ ( m ),
    level_ ( level ),
    index_ ( index ),
    flip_ ( flip )
{
    end_status_ = UuidIntMap();
    for (auto & end : motif_->ends()) {
        end_status_[end->uuid()] = 1;
    }
    connections_ = MotifTreeConnections();
    
}


BasepairOPs
MotifTreeNode::available_ends() {
    BasepairOPs ends;
    int status;
    for( auto const & end : motif_->ends() ) {
        status = end_status_[end->uuid()];
        if (status == 1) { ends.push_back(end); }
    }
    return ends;
}

void
MotifTreeNode::set_end_status(BasepairOP const & end, int status) {
    end_status_[end->uuid()] = status;
}

int
MotifTreeNode::get_end_status(BasepairOP const & end) {
    /*if(end_status_.find(end->uuid()) == end_status_.end()) {
     throw "cannot find end in get_end_status";
     }*/
    return end_status_[end->uuid()];
}

void
MotifTreeNode::add_connection(MotifTreeConnection const & mtc) {
    connections_.push_back(mtc);
}

