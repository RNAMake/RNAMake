//
//  motif_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iomanip>      // std::setw

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


//add functions ////////////////////////////////////////////////////////////////////////////////////



TreeNodeOP<MotifOP>
MotifTree::_get_parent(
    int parent_index) {
    
    auto parent = tree_.last_node();
    
    //catch non existant parent
    try {
        if(parent_index != -1) { parent = tree_.get_node(parent_index); }
    }
    catch(TreeException e) {
        throw MotifTreeException(
            "could not add motif with parent index: " + std::to_string(parent_index) +
            "there is no node with that index");
    }
    
    return parent;
}

Ints
MotifTree::_get_available_parent_end_pos(
    TreeNodeOP<MotifOP> const & parent,
    int parent_end_index) {
    
    auto avail_pos = Ints();
    
    if(parent_end_index != -1) {
        int avail = parent->available_pos(parent_end_index);
        if(!avail) {
            throw MotifTreeException(
                "cannot add motif to tree as the end with index: " +
                std::to_string(parent_end_index) +" you are trying to "
                "add it to is already filled or does not exist");
        }
        
        if(parent_end_index == parent->data()->block_end_add()) {
            throw MotifTreeException(
                "cannot add motif: to tree as the parent_end_index"
                " supplied is blocked see class Motif");
        }
        
        auto name = parent->data()->ends()[parent_end_index]->name();
        if(connections_.in_connection(parent->index(), name)) {
            throw MotifTreeException(
                "cannot add motif to tree as the end "
                "you are trying to add it to is in a connection");
            }
        
            
        avail_pos.push_back(parent_end_index);
    }
    
    else {
        auto avail_pos_temp = parent->available_children_pos();
        for(auto const & p : avail_pos_temp) {
            if(p == parent->data()->block_end_add()) { continue; }
            auto parent_end_name = parent->data()->end_name(p);
            if(connections_.in_connection(parent->index(), parent_end_name)) {
                continue;
            }
            avail_pos.push_back(p);
        }
    }
    
    return avail_pos;
    
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

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    String parent_end_name) {
    
    auto parent = _get_parent(parent_index);
    auto parent_end_index = -1;
    
    try{
        parent_end_index = parent->data()->end_index(parent_end_name);
    }
    catch(RNAStructureException) {
        throw MotifTreeException(
            "cannot find parent_end_name: " + parent_end_name +
            " cannot add motif: " + m->name());
    }
    
    return add_motif(m, parent_index, parent_end_index);
}

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    for(auto const & n : tree_) {
        if(n->data()->id() == m->id()) {
            throw MotifTreeException(
                "cannot add motif there is already a motif in this tree with its unique "
                " indentifier");
        }
    }
    
    auto parent = _get_parent(parent_index);
    
    if(parent == nullptr) {
        auto m_copy = std::make_shared<Motif>(*m);
        m_copy->get_beads(m_copy->ends()[0]);
        int pos = tree_.add_data(m_copy, (int)m_copy->ends().size(), -1, -1);
        merger_.add_motif(m_copy);
        return pos;
    }
    
    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        int pos = tree_.add_data(m_added, (int)m_added->ends().size(), parent->index(), p);
        if(pos != -1) {
            merger_.add_motif(m_added, m_added->ends()[0],
                              parent->data(), parent->data()->ends()[p]);
        }
        return pos;
    }
    
    
    return -1;
    
}

int
MotifTree::_get_connection_end(
    TreeNodeOP<MotifOP> const & node,
    String const & bp_name) {
    
    int node_end_index = -1;
    
    if(bp_name != "") {
        auto ei = node->data()->end_index(bp_name);
        
        if(!node->available_pos(ei)) {
            throw MotifTreeException(
                "cannot add connection with " + std::to_string(node->index()) + " and "
                "end name " + bp_name + " as the end is blocked");
        }
        
        if(ei == node->data()->block_end_add()) {
            throw MotifTreeException(
                "cannot add connection with " + std::to_string(node->index()) + " and "
                "end name " + bp_name + " as the end is blocked");
        }
            
        if(connections_.in_connection(node->index(), bp_name)) {
            throw MotifTreeException(
                "cannot add connection with " + std::to_string(node->index()) +
                " and end name " + bp_name + " as this end is "
                "already in a connection");
        }
        node_end_index = ei;

    }
    
    else {
        auto node_indexes = node->available_children_pos();
        node_indexes.erase(node_indexes.begin(),node_indexes.begin()+1);

        if(node_indexes.size() > 1) {
            throw MotifTreeException(
                "cannot connect nodes " + std::to_string(node->index()) + " its unclear "
                " which ends to attach as there is more then one possibility");
        }
        
        if(node_indexes.size() == 0) {
            throw MotifTreeException(
                "cannot connect nodes " + std::to_string(node->index())  + " there are "
                "no ends free ends to attach too");
        }
        
        auto node_index_name = node->data()->end_name(node_indexes[0]);
        if(connections_.in_connection(node->index(), node_index_name)) {
            throw MotifTreeException(
                "cannot connect nodes " + std::to_string(node->index())  + " there are "
                "no ends free ends to attach too");
        }
        
        node_end_index = node_indexes[0];

    }
    
    return node_end_index;
    
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
    
    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);
    
    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);
    
    
    connections_.add_connection(i, j, node_i_end_name, node_j_end_name);
    
    merger_.connect_motifs(node_i->data(), node_j->data(),
                           node_i->data()->ends()[node_i_ei],
                           node_j->data()->ends()[node_j_ei]);
    
}


//removal functions ////////////////////////////////////////////////////////////////////////////////


void
MotifTree::remove_node(
    int i) {
    
    if(i == -1) {
        i = last_node()->index();
    }
    
    try {
        auto n = get_node(i);
        tree_.remove_node(n);
        merger_.remove_motif(n->data());
        connections_.remove_connections_to(i);
        
    }
    catch(MotifTreeException) {
        throw MotifTreeException(
                                 "cannot remove node with index: " + std::to_string(i) + " as it does not exist");
    }
}

void
MotifTree::remove_node_level(
    int level) {
    
    if(level == -1) { level = tree_.level(); }
    
    auto remove = TreeNodeOPs<MotifOP>();
    for(auto const & n : tree_) {
        if(n->level() >= level) {
            remove.push_back(n);
        }
    }
    
    std::reverse(remove.begin(), remove.end());
    
    for(auto const & n : remove) {
        remove_node(n->index());
    }
    
}


//outputting functions /////////////////////////////////////////////////////////////////////////////


void
MotifTree::write_pdbs(String const & fname) {
    for( auto const & n : tree_) {
        n->data()->to_pdb(fname + "." + std::to_string(n->index()) + ".pdb");
    }
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



//private option functions /////////////////////////////////////////////////////////////////////////


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















