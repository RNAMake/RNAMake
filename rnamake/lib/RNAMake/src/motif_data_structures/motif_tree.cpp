//
//  motif_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_tree.h"
#include "resources/resource_manager.h"

MotifTree::MotifTree(
    String const & s):
tree_(TreeStatic<MotifOP>()),
merger_(MotifMerger()),
options_(Options("MotifTreeOptions"))  {
    setup_options();
    
    auto spl = split_str_by_delimiter(s, "|");
    auto node_spl = split_str_by_delimiter(spl[0], " ");
    int i = -1;
    for(auto const & e : node_spl) {
        i++;
        auto n_spl = split_str_by_delimiter(e, ",");
        auto m = ResourceManager::getInstance().get_motif(n_spl[0], n_spl[2], n_spl[1]);
        if(i == 0) {
            add_motif(m);
        }
        else {
            add_motif(m, std::stoi(n_spl[3]), std::stoi(n_spl[4]));
        }
        
    }
    
    if(spl.size() == 1) { return; }
    
    auto connection_spl = split_str_by_delimiter(spl[1], " ");
    for(auto const & c_str : connection_spl) {
        auto c_spl = split_str_by_delimiter(c_str, ",");
        auto mc = std::make_shared<MotifConnection>(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                                                    c_spl[2], c_spl[3]);
        connections_.push_back(mc);
    }
}

void
MotifTree::setup_options() {
    options_.add_option("sterics", true, OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifTree::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
}

int
MotifTree::add_motif(
    String const & m_name,
    int parent_index,
    String const & p_end_name) {
 
    
    auto parent = last_node();
    if(parent_index != -1) {
        parent = get_node(parent_index);
    }
    auto parent_end_index = parent->data()->end_index(p_end_name);
    return add_motif(m_name, parent_index, parent_end_index);
    
}

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    String parent_end_name) {
    
    auto parent = last_node();
    if(parent_index != -1) {
        parent = get_node(parent_index);
    }
    auto parent_end_index = parent->data()->end_index(parent_end_name);
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifTree::add_motif(
    String const & m_name,
    int parent_index,
    int parent_end_index) {
    
    auto m = MotifOP();
    try {
        m = ResourceManager::getInstance().get_motif(m_name);
    }
    catch(ResourceManagerException const & e) {
        throw MotifTreeException("failed to retrieve motif by name in add_motif: "
                                   + String(e.what()));
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifTree::add_motif(
    String const & m_name,
    String const & m_end_name,
    int parent_index,
    int parent_end_index) {
    
    auto m = MotifOP();
    try {
        m = ResourceManager::getInstance().get_motif(m_name, "", m_end_name);
    }
    catch(ResourceManagerException const & e) {
        throw MotifTreeException("failed to retrieve motif by name in add_motif: "
                                 + String(e.what()));
    }
    
    return add_motif(m, parent_index, parent_end_index);
}


int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    auto parent = tree_.last_node();
    
    //catch out of bounds node index
    try {
        if(parent_index != -1) {
            parent = tree_.get_node(parent_index);
        }
    }
    catch(TreeException e) {
        throw MotifTreeException("could not add motif: " + m->name() + " with parent: "
                                  + std::to_string(parent_index) + "there is no node with" +
                                  "that index");
    }
    
    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->get_beads(m_copy->ends());
        m_copy->new_res_uuids();
        int pos = tree_.add_data(m_copy, (int)m_copy->ends().size(), -1, -1);
        merger_.add_motif(m_copy);
        return pos;
    }
    
    Ints avail_pos;
    try {
        avail_pos = tree_.get_available_pos(parent, parent_end_index);
    }
    catch(TreeException e) {
        throw MotifTreeException("could not add motif: " + m->name() + " with parent: "
                                  + std::to_string(parent_index));
    }
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        m_added->new_res_uuids();
        int pos = tree_.add_data(m_added, (int)m_added->ends().size(), parent->index(), p);
        if(pos != -1) {
            merger_.add_motif(m_added, m_added->ends()[0], parent->data(), parent->data()->ends()[p]);
        }
        return pos;
    }


    return -1;
    
}


void
MotifTree::add_connection(
    int i,
    int j,
    String const & i_bp_name,
    String const & j_bp_name) {
    
    auto node_i = tree_.get_node(i);
    auto node_j = tree_.get_node(j);
    auto name_i = String("");
    auto name_j = String("");
    auto ei = -1;
    auto ej = -1;
    
    if (i_bp_name != "") {
        ei = node_i->data()->end_index(i_bp_name);
        assert(node_i->available_pos(ei) && "cannot add_connection");
        name_i = i_bp_name;
    }
    else {
        auto node_i_indexes = node_i->available_children_pos();
        assert(node_i_indexes.size() != 1 && "cannot add_connection no available spots");
        ei = node_i_indexes[1];
        name_i = node_i->data()->ends()[ei]->name();
    }
    
    
    if (j_bp_name != "") {
        ej = node_j->data()->end_index(j_bp_name);
        assert(node_j->available_pos(ej) && "cannot add_connection");
        name_j = j_bp_name;
    }
    else {
        auto node_j_indexes = node_j->available_children_pos();
        assert(node_j_indexes.size() != 1 && "cannot add_connection no available spots");
        ej = node_j_indexes[1];
        name_j = node_j->data()->ends()[ej]->name();
    }
    
    
    auto connection = std::make_shared<MotifConnection>(i, j, name_i, name_j);
    connections_.push_back(connection);
    
    merger_.connect_motifs(node_i->data(), node_j->data(),
                           node_i->data()->ends()[ei],
                           node_j->data()->ends()[ej]);
    
}


void
MotifTree::write_pdbs(String const & fname) {
    std::stringstream ss;
    for( auto const & n : tree_) {
        ss << fname << "." << n->index() << ".pdb";
        n->data()->to_pdb(ss.str());
        ss.str("");
    }
}

int
MotifTree::_steric_clash(MotifOP const & m) {
    float dist = 0;
    for(auto const & n : tree_) {
        for(auto const & c1 : n->data()->beads()) {
            if(c1.btype() == BeadType::PHOS) { continue; }
            for(auto const & c2 : m->beads()) {
                if(c2.btype() == BeadType::PHOS) { continue; }
                dist = c1.center().distance(c2.center());
                if(dist < clash_radius_) { return 1; }
            }
        }
    }
    return 0;
}

String
MotifTree::topology_to_str() {
    String s;
    
    for(auto const & n : tree_) {
        s += n->data()->name() + "," + n->data()->ends()[0]->name() + ",";
        s += n->data()->end_ids()[0] + "," + std::to_string(n->parent_index()) + ",";
        s += std::to_string(n->parent_end_index()) + " ";
    }
    s += "|";
    for(auto const & c : connections_) {
        s += c->to_str() + " ";
    }
    
    return s;
}


















