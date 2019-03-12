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
#include "structure/residue_type_set_manager.h"


MotifTree::MotifTree():
    tree_(data_structure::tree::TreeStatic<MotifOP>()),
    merger_(nullptr),
    update_merger_(1),
    connections_(MotifConnections()),
    options_(base::Options()) { setup_options(); }

MotifTree::MotifTree(
    String const & s):
    tree_(data_structure::tree::TreeStatic<MotifOP>()),
    merger_(nullptr),
    update_merger_(1),
    connections_(MotifConnections()),
    options_(base::Options())  {
    setup_options();
    
    set_option_value("sterics", false);
    
    auto spl = base::split_str_by_delimiter(s, "|");
    auto node_spl = base::split_str_by_delimiter(spl[0], " ");
    int i = -1;
    int pos = 0;
    for(auto const & e : node_spl) {
        i++;
        auto n_spl = base::split_str_by_delimiter(e, ",");
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
    
    auto connection_spl = base::split_str_by_delimiter(spl[1], " ");
    for(auto const & c_str : connection_spl) {
        auto c_spl = base::split_str_by_delimiter(c_str, ",");
        connections_.add_connection(std::stoi(c_spl[0]), std::stoi(c_spl[1]),
                                    c_spl[2], c_spl[3]);
    }
    
    set_option_value("sterics", true);

}

MotifTree::MotifTree(
        String const & s,
        MotifTreeStringType type):
        tree_(data_structure::tree::TreeStatic<MotifOP>()),
        merger_(nullptr),
        update_merger_(1),
        connections_(MotifConnections()),
        options_(base::Options()) {

    setup_options();
    if(type == MotifTreeStringType::MT_STR) {
        _setup_from_str(s);
    }

}

MotifTree::MotifTree(
    MotifTree const & mt):
    tree_(data_structure::tree::TreeStatic<MotifOP>(mt.tree_)) {
    auto motifs = MotifOPs();
    // dear god this is horrible but cant figure out a better way to do a copy
    for(auto const & n : mt) {
        tree_.get_node(n->index())->data() = std::make_shared<Motif>(*n->data());
        motifs.push_back(tree_.get_node(n->index())->data());
    }
    options_ = base::Options(mt.options_);
    connections_ = MotifConnections(mt.connections_);
    update_merger_ = 1;
}

void
MotifTree::_setup_from_str(
        String const & mt_str) {
    set_option_value("sterics", false);
    auto spl = base::split_str_by_delimiter(mt_str, "FAF");
    auto node_spl = base::split_str_by_delimiter(spl[0], "KAK");
    auto max_index = -100;
    for(auto const & n_str : node_spl) {
        if(n_str.length() < 10) { break; }
        auto n_spl = base::split_str_by_delimiter(n_str, "^");
        auto m = std::make_shared<Motif>(n_spl[0],
                                         structure::ResidueTypeSetManager::getInstance().residue_type_set());
        try {
            auto m2 = RM::instance().motif(m->name());
        } catch(ResourceManagerException const & e) {
            RM::instance().register_motif(m);
        }

        if(m->ends().size() > 1) {
            m->get_beads(m->ends()[0]);
        }
        auto n_index = std::stoi(n_spl[1]);
        tree_.add_data(m, (int)m->ends().size(), std::stoi(n_spl[2]), std::stoi(n_spl[3]));
        tree_.last_node()->index(n_index);

        if(n_index > max_index) {
            max_index = n_index;
        }
    }
    tree_.index(max_index+1);


}


//add functions ////////////////////////////////////////////////////////////////////////////////////

data_structure::tree::TreeNodeOP<MotifOP>
MotifTree::_get_parent(
    int parent_index) {
    
    auto parent = tree_.last_node();
    
    //catch non existant parent
    try {
        if(parent_index != -1) { parent = tree_.get_node(parent_index); }
    }
    catch(data_structure::tree::TreeException const & e) {
        throw MotifTreeException(
            "could not add motif with parent index: " + std::to_string(parent_index) +
            "there is no node with that index");
    }
    
    return parent;
}

Ints
MotifTree::_get_available_parent_end_pos(
    data_structure::tree::TreeNodeOP<MotifOP> const & parent,
    int parent_end_index) {
    
    auto avail_pos = Ints();
    
    if(parent == nullptr) { return avail_pos; }
    
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
            if(c1.btype() == structure::BeadType::PHOS) { continue; }
            for(auto const & c2 : m->beads()) {
                if(c2.btype() == structure::BeadType::PHOS) { continue; }
                dist = c1.center().distance(c2.center());
                if(dist < clash_radius_) { return 1; }
            }
        }
    }
    return 0;
}

int
MotifTree::_get_connection_end(
    data_structure::tree::TreeNodeOP<MotifOP> const & node,
    String const & bp_name) {
    
    int node_end_index = -1;
    
    if(bp_name != "") {
        auto ei = node->data()->get_end_index(bp_name);
        
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

int
MotifTree::add_motif(
    MotifOP const & m,
    int parent_index,
    String parent_end_name) {
    
    auto parent = _get_parent(parent_index);
    auto parent_end_index = -1;
    
    try{
        parent_end_index = parent->data()->get_end_index(parent_end_name);
    }
    catch(structure::RNAStructureException) {
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
        return pos;
    }
    
    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        auto m_added = get_aligned_motif(parent->data()->ends()[p], m->ends()[0], m);
        if(sterics_ && _steric_clash(m_added)) { continue; }
        
        int pos = tree_.add_data(m_added, (int)m_added->ends().size(), parent->index(), p);
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
    
    auto node_i = data_structure::tree::TreeNodeOP<MotifOP>(nullptr);
    auto node_j = data_structure::tree::TreeNodeOP<MotifOP>(nullptr);

    try {  node_i = tree_.get_node(i); }
    catch(data_structure::tree::TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(i) +" does not exist");
    }
    
    try {  node_j = tree_.get_node(j); }
    catch(data_structure::tree::TreeException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(j) +" does not exist");
    }
    
    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);
    
    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);
    
    
    connections_.add_connection(i, j, node_i_end_name, node_j_end_name);
    
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
    
    auto remove = data_structure::tree::TreeNodeOPs<MotifOP>();
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

String
MotifTree::to_str() {
    auto s = String();
    for(auto const & n : tree_) {
        s += n->data()->to_str() + "^" + std::to_string(n->index()) + "^";
        s += std::to_string(n->parent_index()) + "^" + std::to_string(n->parent_end_index()) + " KAK ";
    }
    s += " FAF ";
    for(auto const & c : connections_) {
        s += c->to_str() + "|";
    }
    return s;
}


//misc /////////////////////////////////////////////////////////////////////////////////////////////

void
MotifTree::_update_merger() {
    if(! update_merger_) { return; }
    
    merger_ = std::make_shared<MotifMerger>();
    
    for(auto const & n : tree_) {
        if(n->parent() == nullptr) {
            merger_->add_motif(n->data());
            continue;
        }
        
        auto pei = n->parent_end_index();
        merger_->add_motif(n->data(), n->data()->ends()[0],
                           n->parent()->data(), n->parent()->data()->ends()[pei]);
    }
    
    for(auto const & c : connections_) {
        auto node_i = get_node(c->i());
        auto node_j = get_node(c->j());
        auto node_i_ei = node_i->data()->get_end_index(c->name_i());
        auto node_j_ei = node_j->data()->get_end_index(c->name_j());
        merger_->connect_motifs(node_i->data(), node_j->data(),
                                node_i->data()->ends()[node_i_ei],
                                node_j->data()->ends()[node_j_ei]);
    }
    
    update_merger_ = 0;
}


//getters //////////////////////////////////////////////////////////////////////////////////////////


structure::Beads
MotifTree::beads() {
    auto beads = structure::Beads();
    for(auto const & n : tree_) {
        for(auto b : n->data()->beads()) {
            beads.push_back(b);
        }
    }
    return beads;
}



//private option functions /////////////////////////////////////////////////////////////////////////


void
MotifTree::setup_options() {
    options_.add_option("sterics", true, base::OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, base::OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifTree::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
}















