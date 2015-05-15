//
//  motif_tree_node.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <algorithm>

//RNAMake Headers
#include "base/types.h"
#include "util/settings.h"
#include "structure/residue_type_set.h"
#include "motif/motif.h"
#include "motif/motif_tree_node.h"


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
        end_status_[end->uuid().s_uuid()] = 1;
    }
    connections_ = MotifTreeConnectionOPs();
    
}

MotifTreeNode::~MotifTreeNode() {
}

BasepairOPs
MotifTreeNode::available_ends() {
    BasepairOPs ends;
    int status;
    for( auto const & end : motif_->ends() ) {
        status = end_status_[end->uuid().s_uuid()];
        if (status == 1) { ends.push_back(end); }
    }
    return ends;
}

void
MotifTreeNode::set_end_status(BasepairOP const & end, int status) {
    end_status_[end->uuid().s_uuid()] = status;
}

int
MotifTreeNode::get_end_status(BasepairOP const & end) {
    /*if(end_status_.find(end->uuid()) == end_status_.end()) {
     throw "cannot find end in get_end_status";
     }*/
    return end_status_[end->uuid().s_uuid()];
}

void
MotifTreeNode::add_connection(MotifTreeConnectionOP const & mtc) {
    connections_.push_back(mtc);
}

void
MotifTreeNode::remove_connection(MotifTreeConnectionOP const & c) {

    try {
        connections_.erase(std::remove(connections_.begin(), connections_.end(), c), connections_.end());
    }
    catch(...) {
        throw "could not remove connection in motif_tree_node"; 
    }
}

MotifTreeNodeOP
MotifTreeNode::parent() {
    MotifTreeNodeOP partner;
    for(auto const & c : connections_) {
        if (c->node_1().get() == this ) {
            partner = c->node_2();
        }
        else {
            partner = c->node_1();
        }
        
        if(partner->index() < index_) { return partner; }
    }
    
    throw "something crazy happened, couldnt find the parent of motif tree node";
}

MotifTreeConnectionOP const &
MotifTreeNode::connection(
    MotifTreeNodeOP const & node) {
    for(auto const &  c : connections_) {
        if(c->node_1() == node || c->node_2() == node) {
            return c;
        }
    }
    throw "connection called with node that is not connected to current node";
}

int
MotifTreeNode::motif_end_index(
    MotifTreeNodeOP const & node) {
    
    MotifTreeConnectionOP c = connection(node);
    BasepairOP end = c->motif_end(node);
    return motif_->end_index(end);
}

int
MotifTreeNode::parent_end_index() {
    return motif_end_index(parent());
}




















