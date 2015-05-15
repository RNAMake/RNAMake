//
//  motif_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <algorithm>

//RNAMake Headers
#include "util/settings.h"
#include "structure/resource_manager.h"
#include "motif/motif_tree.h"
#include "motif/motif.h"

MotifTree::MotifTree() {
    String path = resources_path() + "/start.motif";
    String line;
    std::ifstream in;
    in.open(path.c_str());
    getline(in, line);
    in.close();
    MotifOP m ( new Motif ( line, ResourceManager::getInstance().residue_type_set()) );
    MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
    last_node_ = head;
    nodes_ = MotifTreeNodeOPs();
    nodes_.push_back(head);
    clash_radius_ = 2.8;
    sterics_ = 1;
    full_beads_first_res_ = 1;
    merger_ = MotifTreeMerger();
    level_ = 1;
}

MotifTree::MotifTree(MotifOP const & m) {
    MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
    last_node_ = head;
    nodes_ = MotifTreeNodeOPs();
    nodes_.push_back(head);
    clash_radius_ = 2.8;
    sterics_ = 1;
    full_beads_first_res_ = 1;
    merger_ = MotifTreeMerger();
    level_ = 1;

}

MotifTree::MotifTree(
    String const & s,
    ResidueTypeSet const & rts) {
    nodes_ = MotifTreeNodeOPs();
    clash_radius_ = 2.8;
    sterics_ = 1;
    full_beads_first_res_ = 1;
    Strings spl = split_str_by_delimiter(s, "#");
    int i = -1;
    for (auto const & e : spl) {
        i++;
        Strings mtn_spl = split_str_by_delimiter(e, "!");
        MotifOP m ( new Motif ( mtn_spl[0], rts));
        if ( i == 0) {
            MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
            nodes_.push_back(head);
            continue;
        }
        MotifTreeNodeOP mtn ( new MotifTreeNode(m, 1, i, std::stoi(mtn_spl[4])));
        MotifTreeNodeOP parent = nodes_[ std::stoi(mtn_spl[1]) ];
        BasepairOP parent_end = parent->motif()->ends()[ std::stoi(mtn_spl[2]) ];
        BasepairOP node_end = m->ends() [ std::stoi(mtn_spl[3]) ];
        MotifTreeConnection(parent, mtn, parent_end, node_end);
        nodes_.push_back(mtn);
    }
    last_node_ = nodes_.back();
    merger_ = MotifTreeMerger();
    level_ = 1;
}


MotifTreeNodeOP
MotifTree::add_motif(
    MotifOP const & m,
    MotifTreeNodeOP parent,
    int end_index,
    int parent_index,
    int end_flip) {
    
    if(parent == nullptr) { parent = last_node_; }
    
    BasepairOPs parent_ends, motif_ends;
    
    if (parent_index == -1) {
        parent_ends = parent->available_ends();
        if(parent_ends.size() == 0) {
            throw "cannot use this parent it has not available ends\n";
        }
    
    }
    else                    {
        if(parent_index >= parent->motif()->ends().size()) {
            throw "cannot use this parent index, it is out of range\n";
        }
        int status = parent->get_end_status(parent->motif()->ends()[parent_index]);
        if(status == 0) {
            throw "cannot ue this parent index, it is currently being used";
        }
        parent_ends.push_back(parent->motif()->ends()[parent_index]);
    }
    
    if (end_index == -1) {
        
        motif_ends = m->ends();
        if(motif_ends.size() == 0) {
            throw "cannot add this motif to the tree it has no ends\n";
        }
    }
    else                 {
        if(end_index >= m->ends().size()) {
            throw "cannot use this end index, it is out of range\n";
        }
        motif_ends.push_back(m->ends()[end_index]);
    }

    MotifTreeNodeOP new_node(NULL);
    Ints flips;
    if ( end_flip == -1) {
        flips.push_back(0);
        flips.push_back(1);
    }
    else {
        if(end_flip != 0 && end_flip != 1) {
            throw "supplied an end_flip that is not 1 or 0\n";
        }
        flips.push_back(end_flip);
    }

    for (auto const & parent_end : parent_ends) {
        for (auto const & end : motif_ends) {
            new_node = _add_motif(parent, m, parent_end, end, flips);
            if( new_node.get() != NULL) { break; }
        }
        if ( new_node.get() != NULL ) { break; }
    }
    
    if ( new_node.get() != NULL ) {
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
        end->flip(flip);
        align_motif(parent_end, end, m);
        if(sterics_ == 1) {
            if(_steric_clash(m, end) == 1) {
                m->reset();
                continue;
            }
        }
        
        MotifOP m_copy ( new Motif ( m->copy() ));
        MotifTreeNodeOP new_node ( new MotifTreeNode(m_copy, level_, (int)nodes_.size(), flip) );
        m->reset();
        BasepairOP new_end = m_copy->get_basepair(end->uuid())[0];
        for ( auto const & r : m_copy->residues() ) {  r->new_uuid(); }
        MotifTreeConnectionOP mtc (new MotifTreeConnection(parent, new_node, parent_end, new_end, 0));
        parent->add_connection(mtc);
        new_node->add_connection(mtc);
        
        if(nodes_.size() == 1 && full_beads_first_res_ == 1) {
            nodes_[0]->motif()->get_beads();
        }
        
        _update_beads(parent, new_node);
        return new_node;
    }
    
    return MotifTreeNodeOP(NULL);
    
}

int
MotifTree::_steric_clash(
    MotifOP const & m,
    BasepairOP const & end) {
    
    Beads beads = m->get_beads(end);
    float dist;
    for( auto const & n : nodes_) {
        for ( auto const & b1 : n->motif()->beads() ) {
            for ( auto const & b2 : beads ) {
                if (b1.btype() == PHOS || b2.btype() == PHOS) { continue; }
                dist = b1.center().distance(b2.center());
                if(dist < clash_radius_) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

void
MotifTree::_update_beads(
    MotifTreeNodeOP const & parent,
    MotifTreeNodeOP const & child) {
    
    if (parent->motif()->mtype() != HELIX || child->motif()->mtype() != HELIX) {
        return;
    }
    
    BasepairOPs not_used = parent->available_ends();
    BasepairOPs exclude;
    int found = 0;
    for( auto const & end1 : parent->motif()->ends()) {
        found = 0;
        for ( auto const & end2 : not_used) {
            if(end1->uuid() == end2->uuid())  { found = 1; break; }
        }
        if(!found) { exclude.push_back(end1); }
    }
    parent->motif()->get_beads(exclude);
    child->motif()->get_beads();
}


void
MotifTree::_add_connection(
    MotifTreeNodeOP const & node_1,
    MotifTreeNodeOP const & node_2,
    float cutoff) {
    
    BasepairOPs avail_ends_1 = node_1->available_ends();
    BasepairOPs avail_ends_2 = node_2->available_ends();
    float dist;

    for(auto const & end1 : avail_ends_1) {
        for(auto const & end2 : avail_ends_2) {
            dist = end1->d().distance(end2->d());
            if(dist < cutoff) {
                MotifTreeConnectionOP mtc (new MotifTreeConnection(node_1, node_2, end1, end2));
                node_1->add_connection(mtc);
                node_2->add_connection(mtc);
            }
        }
    }
    
}


void
MotifTree::remove_node(
    MotifTreeNodeOP const & node) {
        
    if(node->index() == 0) {
        throw "cannot remove head node from tree";
    }
    
    MotifTreeNodeOP parent = node->parent();
    
    for(auto & c : node->connections()) {
        MotifTreeNodeOP partner = c->partner(node);
        partner->remove_connection(c);
        c->disconnect();
    }
    
    try {
        nodes_.erase(std::remove(nodes_.begin(), nodes_.end(), node), nodes_.end());
    }
    catch(...) {
        throw "could not remove motif_tree_node from motif_tree";
    }
    
    last_node_ = parent;
    if(nodes_.size() == 1) {
        nodes_[0]->motif()->remove_beads();
    }
}

void
MotifTree::write_pdbs(String const & fname) {
    int i = 0;
    std::stringstream ss;
    for( auto const & n : nodes_) {
        ss << fname << "." << i << ".pdb";
        n->motif()->to_pdb(ss.str());
        ss.str("");
        i++;
    }
}

PoseOP
MotifTree::to_pose(int include_head) {
    return merger_.merge(*this);
}

void
MotifTree::to_pdb(String fname, int include_head) {
    PoseOP pose = to_pose();
    pose->to_pdb(fname);
}


