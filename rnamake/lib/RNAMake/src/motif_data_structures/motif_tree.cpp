//
//  motif_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_tree.h"
#include "resources/resource_manager.h"

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
        merger_.add_motif(m_added, m_added->ends()[0], parent->data(), parent->data()->ends()[p]);
        return pos;
    }


    return -1;
    
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
