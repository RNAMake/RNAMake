//
//  motif_state_ensemble_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "data_structure/graph/graph_node.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"
#include "motif/motif_ensemble.h"
#include "resources/resource_manager.h"


int
MotifStateEnsembleTree::add_ensemble(
    MotifStateEnsembleOP const & ensemble,
    int parent_index,
    int parent_end_index) {
    
    auto parent = tree_.last_node();
    if(parent_index != -1) {
        parent = tree_.get_node(parent_index);
    }
    
    if(parent == nullptr) {
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(ensemble->copy()),
                              ensemble->num_end_states());
    }
    
    auto avail_pos = tree_.get_available_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        if(p == parent->data()->block_end_add()) { continue; }
        
        return tree_.add_data(std::make_shared<MotifStateEnsemble>(ensemble->copy()),
                              ensemble->num_end_states(),
                              parent->index(),
                              p);
    }
    
    return -1;
    
}

MotifStateTreeOP
MotifStateEnsembleTree::to_mst() {
    
    auto mst = std::make_shared<MotifStateTree>();
    mst->option("sterics", 0);
    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : tree_) {
        i++;
        auto state = n->data()->most_populated();
        if(i == 0) {
            mst->add_state(state);
            continue;
        }
        
        parent_index = n->parent_index();
        parent_end_index = n->parent_end_index();
        j = mst->add_state(state, parent_index, parent_end_index);
        if(j == -1) {
            std::runtime_error("can not build motif state tree from mset");
        }
    }
    
    return mst;
}

void
MotifStateEnsembleTree::setup_from_mt(
    MotifTreeOP const & mt) {
    
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : *mt) {
        i++;
        MotifStateEnsembleOP mse;
    
        if(n->data()->mtype() == HELIX) {
            mse = ResourceManager::getInstance().get_motif_state_ensemble(n->data()->end_ids()[0]);
        }
        
        else {
            int found = ResourceManager::getInstance().has_supplied_motif_state_ensemble(
                            n->data()->name(),
                            n->data()->ends()[0]->name());
            if(found) {
                mse = ResourceManager::getInstance().get_registered_extra_motif_state_ensemble(
                                                            n->data()->name(),
                                                            n->data()->ends()[0]->name());
            }
            
            else {
                auto m = ResourceManager::getInstance().get_motif(n->data()->name(),
                                                                  n->data()->end_ids()[0]);
                mse = std::make_shared<MotifStateEnsemble>(m->get_state());
            }
            
        }
        
        if(i == 0) {
            add_ensemble(mse);
        }
        else {
            parent_index = n->parent()->index();
            parent_end_index = n->parent_end_index();
            if(parent_end_index == -1) {
                std::runtime_error("cannot setup_from_mt in MotifStateEnsembleTree");
            }
            j = add_ensemble(mse, parent_index, parent_end_index);
            if(j == -1) {
                std::runtime_error("failed to add ensemble in setup_from_mt");
            }
        }
        
        
    }
    
    
}



























