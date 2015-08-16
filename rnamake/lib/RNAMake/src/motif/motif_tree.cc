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
#include "motif/motif_tree.h"
#include "motif/motif.h"



MotifTree::MotifTree() {
    graph_ = GraphStatic<MotifOP>();
    /*MotifOP m ( new Motif ( line, ResidueTypeSetManager::getInstance().residue_type_set()) );
    MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
    last_node_ = head;
    nodes_ = MotifTreeNodeOPs();
    nodes_.push_back(head);
    merger_ = MotifTreeMerger();
    level_ = 1;*/
    setup_options(); update_var_options();
}

MotifTree::MotifTree(MotifOP const & m) {
    /*MotifTreeNodeOP head ( new MotifTreeNode (m, 0, 0, 0));
    last_node_ = head;
    nodes_ = MotifTreeNodeOPs();
    nodes_.push_back(head);
    merger_ = MotifTreeMerger();
    level_ = 1;
    setup_options(); update_var_options();*/

}


void
MotifTree::setup_options() {
    options_ = Options();
    options_.add_option(Option("sterics", 1));
    options_.add_option(Option("clash_radius", 2.9f));
}

void
MotifTree::update_var_options() {
    sterics_              = options_.option<int>("sterics");
    clash_radius_         = options_.option<float>("clash_radius");
}

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index,
    String parent_end_name) {
    
    auto parent = graph_.last_node();
    
    //catch out of bounds node index
    try {
        if(parent_index != -1) { parent = graph_.get_node(parent_index); }
    }
    catch(GraphException e) {
        throw MotifTreeException("could not add motif: " + m->name() + " with parent: " + std::to_string(parent_index) + "there is no node with that index");
    }
    
    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(m->copy());
        m_copy->get_beads(m_copy->ends());
        return graph_.add_data(m_copy, -1, -1, -1, (int)m_copy->ends().size());
        
    }
    
    if(parent_end_name.size() > 0) {
        auto parent_end = parent->data()->get_basepair(parent_end_name)[0];
        parent_end_index = parent->data()->end_index(parent_end);
    }
    
    //catch not a viable end index for this parent node
    Ints avail_pos;
    try {
        avail_pos = graph_.get_available_pos(parent, parent_end_index);
    }
    catch(GraphException e) {
        throw MotifTreeException("could not add motif: " + m->name() + " with parent: " + std::to_string(parent_index));
    }
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        m_added->new_res_uuids();
        return graph_.add_data(m_added, parent->index(), p, 0, (int)m_added->ends().size());

     
    }
    
    return -1;
}

int
MotifTree::_steric_clash(
    MotifOP const & m) {
    
    float dist;
    for( auto const & n : graph_) {
        for ( auto const & b1 : n->data()->beads() ) {
            for ( auto const & b2 : m->beads() ) {
                if (b1.btype() == BeadType::PHOS || b2.btype() == BeadType::PHOS) { continue; }
                dist = b1.center().distance(b2.center());
                if(dist < clash_radius_) {
                    return 1;
                }
            }
        }
    }
    return 0;
}



/*void
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
    
}*/


void
MotifTree::remove_node(
    int i) {

    if(graph_.size() == 0) {
        throw MotifTreeException("tried to remove a node from motiftree but there arent any nodes!");
    }
    
    try {
        if(i == -1) {
            i = graph_.last_node()->index();
        }
        graph_.remove_node(i);
    }
    catch(GraphException e) {
        throw MotifTreeException("tried to remove node from motiftree that does not exist");
    }
}


void
MotifTree::remove_node_level(int level) {
    
    if(level == -1) { level = graph_.level(); }
    int size = (int)graph_.size();
    
    Ints indices;
    for(auto const & n : graph_) {
        if(n->level() >= level) { indices.push_back(n->index()); }
    }
    
    for(auto & i: indices) { graph_.remove_node(i); }
    
}


void
MotifTree::write_pdbs(String const & fname) {
    int i = 0;
    std::stringstream ss;
    for( auto const & n : graph_) {
        ss << fname << "." << i << ".pdb";
        n->data()->to_pdb(ss.str());
        ss.str("");
        i++;
    }
}

/*PoseOP
MotifTree::to_pose(int include_head) {
    return merger_.merge(*this);
}*/

void
MotifTree::to_pdb(String fname, int include_head) {
    /*PoseOP pose = to_pose();
    pose->to_pdb(fname);*/
}




















