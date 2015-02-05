//
//  motif_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_tree.h"
#include "motif.h"
#include "types.h"
#include "settings.h"
#include "residue_type_set.h"

MotifTree::MotifTree() {
    String path = resources_path() + "/start.motif";
    String line;
    std::ifstream in;
    in.open(path.c_str());
    getline(in, line);
    in.close();
    ResidueTypeSet rts;
    MotifOP m ( new Motif ( line, rts) );
    MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
    last_node_ = head;
    nodes_ = MotifTreeNodeOPs();
    nodes_.push_back(head);
}

MotifTree::MotifTree(MotifOP m) {
    
}

MotifTreeNodeOP
MotifTree::add_motif(
    MotifOP const & m,
    MotifTreeNodeOP parent = NULL,
    int end_index = -1,
    int parent_index = -1) {
    
    if(parent == NULL) { parent = last_node_; }
    
    BasepairOPs parent_ends, motif_ends;
    
    if (parent_index == -1) {  parent_ends = parent->available_ends(); }
    else                    {
        if(parent_index >= parent->motif()->ends().size()) {
            throw "cannot use this parent index, it is out of range";
        }
        int status = parent->get_end_status(parent->motif()->ends()[parent_index]);
        if(status == 0) {
            throw "cannot ue this parent index, it is currently being used";
        }
        parent_ends.push_back(parent->motif()->ends()[parent_index]);
    }
    
    if (end_index == -1) {  motif_ends = m->ends(); }
    else                 {
        if(end_index >= m->ends().size()) {
            throw "cannot use this end index, it is out of range";
        }
        motif_ends.push_back(m->ends()[end_index]);
    }

    MotifTreeNodeOP new_node = NULL;
    Ints flips(1);
    flips[0] = 0;

    for (auto const & parent_end : parent_ends) {
        for (auto const & end : motif_ends) {
            new_node = _add_motif(parent, m, parent_end, end, flips);
            if( new_node != NULL) { break; }
        }
        if ( new_node != NULL ) { break; }
    }
    
    if ( new_node != NULL ) {
        nodes_.push_back(new_node);
        last_node_ = new_node;
    }
    
    return new_node;
    
}

MotifTreeNodeOP
MotifTree::_add_motif(
    MotifTreeNodeOP const & parent,
    MotifOP const & m,
    BasepairOP const & parent_end,
    BasepairOP const & end,
    Ints const & flip_status) {
    
    for (auto const & flip : flip_status) {
        align_motif(parent_end, end, m);
        MotifOP m_copy ( new Motif ( m->copy() ));
        MotifTreeNodeOP new_node ( new MotifTreeNode(m_copy, 1, (int)nodes_.size(), flip) );
        m->reset();
        BasepairOP new_end = m_copy->get_basepair(end->uuid())[0];
        for ( auto const & r : m_copy->residues() ) {  r->new_uuid(); }
        MotifTreeConnection(parent, new_node, parent_end, new_end, 0);
        return new_node;
    }
    
    return NULL;
    
}

void
MotifTree::write_pdbs() {
    int i = 0;
    std::stringstream ss;
    for( auto const & n : nodes_) {
        ss << "nodes." << i << ".pdb";
        n->motif()->to_pdb(ss.str());
        ss.str("");
        i++;
    }
}

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


