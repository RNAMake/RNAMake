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
    merger_ = MotifTreeMerger();
    level_ = 0;
    setup_options(); update_var_options();
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


void
MotifTree::add_connection(
    int i,
    int j,
    String const & i_name) {
    
    auto node_1 = get_node(i);
    auto node_2 = get_node(j);
    if (i == j) {
        throw MotifTreeException("cannot connect motif to itself in MotifTree::add_connection");
    }
    
    auto avail_ends_1 = node_1->available_children_pos();
    auto avail_ends_2 = node_2->available_children_pos();

    if(avail_ends_1.size() == 0) {
        throw MotifTreeException("cannot connect motifs together in MotifTree::add_connection node: " + std::to_string(i) + " has no available connections");
    }
    
    if(avail_ends_2.size() == 0) {
        throw MotifTreeException("cannot connect motifs together in MotifTree::add_connection node: " + std::to_string(j) + " has no available connections");
    }
    
    for(auto const & end_index_1 : avail_ends_1) {
        if(i_name != node_1->data()->ends()[end_index_1]->name()) {
            continue;
        }
        
        graph_.connect(i, j, end_index_1, avail_ends_2[0]);
        return;
    }
    
    throw MotifTreeException("could not connect node " + std::to_string(i) + " " + std::to_string(j) + " with end name " + i_name);
    
}


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
    
    Ints indices;
    for(auto const & n : graph_) {
        if(n->level() >= level) { indices.push_back(n->index()); }
    }
    
    for(auto & i: indices) { graph_.remove_node(i); }
    
}

String
MotifTree::topology_to_str() {
    String s = "";
    for(auto const & n : graph_) {
        s += n->data()->name() + "," + n->data()->end_ids()[0] + ",";
        //s += std::to_string(n->parent_index());
    }
    
    return s;
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

PoseOP
MotifTree::to_pose() {
    return merger_.merge(graph_);
}

void
MotifTree::to_pdb(String fname) {
    PoseOP pose = to_pose();
    pose->to_pdb(fname);
}




















