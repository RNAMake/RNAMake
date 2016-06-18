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
connections_(MotifConnections()),
options_(Options())  {
    setup_options();
    
    set_option_value("sterics", false);
    
    auto spl = split_str_by_delimiter(s, "|");
    auto node_spl = split_str_by_delimiter(spl[0], " ");
    int i = -1;
    int pos = 0;
    for(auto const & e : node_spl) {
        i++;
        auto n_spl = split_str_by_delimiter(e, ",");
        auto m = MotifOP(nullptr);
        
        try {
            m = RM::instance().motif(n_spl[0], n_spl[2], n_spl[1]);
        }
        catch(ResourceManagerException const & e) {
            throw MotifTreeException(
                String("could not get motif did you forget to add it to the resource manager: ") +
                e.what());
        }
            
        if(i == 0) {
            pos = add_motif(m);
        }
        else {
            pos = add_motif(m, std::stoi(n_spl[3]), std::stoi(n_spl[4]));
        }
        
        if(pos == -1) {
            throw MotifTreeException(
                "failed to add " + m->name() + " pos " + std::to_string(i) + " in the tree "
                "during rebuild from string");
        }
        
    }
    
    if(spl.size() == 1) { return; }
    
    auto connection_spl = split_str_by_delimiter(spl[1], " ");
    for(auto const & c_str : connection_spl) {
        auto c_spl = split_str_by_delimiter(c_str, ",");
        connections_.add_connection(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                                    c_spl[2], c_spl[3]);
    }
    
    set_option_value("sterics", true);

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
    MotifOP const & m,
    int parent_index,
    String parent_end_name) {
    
    auto parent = _get_parent(m->name(), parent_index);
    auto parent_end_index = -1;
    
    try{
        parent_end_index = parent->data()->end_index(parent_end_name);
    }
    catch(RNAStructureException) {
        throw MotifTreeException(
            "cannot find parent_end_name: " + parent_end_name + " cannot add motif: " + m->name());
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    auto parent = _get_parent(m->name(), parent_index);
    
    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->get_beads(m_copy->ends());
        m_copy->new_res_uuids();
        int pos = tree_.add_data(m_copy, (int)m_copy->ends().size(), -1, -1);
        merger_.add_motif(m_copy);
        return pos;
    }
    
    auto avail_pos = _get_available_parent_end_pos(m->name(), parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        if(connections_.in_connection(parent->index(),
                                      parent->data()->ends()[p]->name())) {
            continue;
        }
        
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
    
    auto node_i = TreeNodeOP<MotifOP>(nullptr);
    auto node_j = TreeNodeOP<MotifOP>(nullptr);

    try {  node_i = tree_.get_node(i); }
    catch(TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(i) +" does not exist");
    }
    
    try {  node_j = tree_.get_node(j); }
    catch(TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(j) +" does not exist");
    }
    
    auto name_i = String("");
    auto name_j = String("");
    auto ei = -1;
    auto ej = -1;
    
    if (i_bp_name != "") {
        ei = node_i->data()->end_index(i_bp_name);
        
        if(!node_i->available_pos(ei)) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(i) + " and " +
                std::to_string(j) + " as specified end for i is filled: " + i_bp_name);
        }
        
        if(connections_.in_connection(i, i_bp_name)) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(i) + " and " +
                std::to_string(j) + " as specified end for is already in connection: " + i_bp_name);
        }
        
        name_i = i_bp_name;
    }
    else {
        auto node_i_indexes = node_i->available_children_pos();
        
        if(node_i_indexes.size() == 1) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(i) + " and " +
                std::to_string(j) + " as no ends are available for : " + std::to_string(i));
        }

        
        ei = node_i_indexes[1];
        name_i = node_i->data()->ends()[ei]->name();
    }
    
    
    if (j_bp_name != "") {
        ej = node_j->data()->end_index(j_bp_name);

        if(!node_j->available_pos(ej)) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(i) + " and " +
                    std::to_string(j) + " as specified end for j is filled: " + j_bp_name);
        }
        
        if(connections_.in_connection(j, j_bp_name)) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(j) + " and " +
                std::to_string(j) + " as specified end for is already in connection: " + j_bp_name);
        }
        name_j = j_bp_name;
    }
    else {
        auto node_j_indexes = node_j->available_children_pos();
        
        if(node_j_indexes.size() == 1) {
            throw MotifTreeException(
                "cannot add connection between nodes " + std::to_string(i) + " and " +
                std::to_string(j) + " as no ends are available for : " + std::to_string(j));
        }
        ej = node_j_indexes[1];
        name_j = node_j->data()->ends()[ej]->name();
    }
    
    
    if(connections_.in_connection(i, node_i->data()->ends()[ei]->name())) {
        throw MotifTreeException(
            "cannot add connection between nodes " + std::to_string(i) + " and " +
            std::to_string(j) + "as only end available for " + std::to_string(i) +
            " is already in a connection already");
    }
    
    if(connections_.in_connection(j, node_j->data()->ends()[ej]->name())) {
        throw MotifTreeException(
            "cannot add connection between nodes " + std::to_string(i) + " and " +
            std::to_string(j) + "as only end available for " + std::to_string(j) +
            " is already in a connection already");
    }
    
    
    connections_.add_connection(i, j, name_i, name_j);
    
    merger_.connect_motifs(node_i->data(), node_j->data(),
                           node_i->data()->ends()[ei],
                           node_j->data()->ends()[ej]);
    
}


void
MotifTree::write_pdbs(String const & fname) {
    for( auto const & n : tree_) {
        n->data()->to_pdb(fname + "." + std::to_string(n->index()) + ".pdb");
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


















