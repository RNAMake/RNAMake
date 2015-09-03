//
//  motif_state_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_state_tree.h"
#include "resources/resource_manager.h"


void
MotifStateTree::setup_options() {
    options_ = Options();
    options_.add_option(Option("sterics", 1));
    options_.add_option(Option("clash_radius", 2.9f));
}

void
MotifStateTree::update_var_options() {
    sterics_              = options_.option<int>("sterics");
    clash_radius_         = options_.option<float>("clash_radius");
}

int
MotifStateTree::add_state(
    MotifStateOP const & state,
    int parent_index,
    int parent_end_index,
    String parent_end_name) {
    
    auto parent = tree_.last_node();

    //catch out of bounds node index
    try {
        if(parent_index != -1) { parent = tree_.get_node(parent_index); }
    }
    catch(TreeException e) {
        throw MotifStateTreeException("could not add state: " + state->name() + " with parent: " + std::to_string(parent_index) + "there is no node with that index");
    }
    
    if(parent == nullptr) {
        auto n_data = std::make_shared<MSTNodeData>(state);
        return tree_.add_data(n_data, (int)state->end_states().size(), -1, -1);
    }
    
    if(parent_end_name.length() > 0) {
        auto parent_end = parent->data()->cur_state->get_end_state(parent_end_name);
        parent_end_index = parent->data()->cur_state->end_index(parent_end);
    }
    
    auto avail_pos = tree_.get_available_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->ref_state->block_end_add()) { continue; }
        auto n_data = std::make_shared<MSTNodeData>(state);
        
        get_aligned_motif_state(parent->data()->cur_state->end_states()[p],
                                n_data->cur_state,
                                n_data->ref_state);
        
        if(sterics_ && _steric_clash(n_data)) { continue; }
    
        return tree_.add_data(n_data, (int)state->end_states().size(), parent->index(), p);
        
    }
    
    return -1;
}

int
MotifStateTree::setup_from_mt(
    MotifTreeOP const & mt) {
    
    int i = -1, j = -1;
    int parent_index = 1000, parent_end_index = -1;
    for(auto const & n : *mt) {
        i++;
        auto ms = ResourceManager::getInstance().get_state(n->data()->name(),
                                                           n->data()->end_ids()[0],
                                                           n->data()->ends()[0]->name());
        if(i == 0) {
            add_state(ms);
        }
        else {
            parent_index = n->parent()->index();
            parent_end_index = n->parent_end_index();
            if(parent_index == -1 || parent_end_index == -1) {
                throw std::runtime_error("did not convert motif tree to motif state tree properly");
            }
            j = add_state(ms, parent_index, parent_end_index);
            if(j == -1) {
                throw std::runtime_error("could not add motif state to tree during conversion");
            }
        }
    }
    
    return 1;
    
}





















