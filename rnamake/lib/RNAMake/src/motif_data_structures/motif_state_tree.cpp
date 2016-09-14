//
//  motif_state_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif_data_structures/motif_state_tree.h"
#include "resources/resource_manager.h"


MotifStateTree::MotifStateTree():
    tree_(TreeStatic<MSTNodeDataOP>()),
    aligner_(MotifStateAligner()),
    queue_(std::queue<MotifStateTreeNodeOP>()),
    options_(Options()) {
    setup_options();
}

MotifStateTree::MotifStateTree(
    MotifTreeOP const & mt):
    MotifStateTree() {
    
    int i = -1;
    for(auto const & n : *mt) {
        i++;
        auto ms = RM::instance().motif_state(n->data()->name(),
                                             n->data()->end_ids()[0],
                                             n->data()->ends()[0]->name());
        
        if(i == 0) {
            add_state(ms);
        }
        
        else {
            int j = add_state(ms, n->parent_index(), n->parent_end_index());
            if(j == -1) {
                throw MotifStateTreeException(
                    "could not convert motif tree to motif state tree");
            }
        }
        
    }
    
}


//add functions ////////////////////////////////////////////////////////////////////////////////////


void
MotifStateTree::setup_options() {
    options_.add_option("sterics", true, OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifStateTree::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
}

int
MotifStateTree::add_state(
    MotifStateOP const & state,
    int parent_index,
    int parent_end_index,
    String parent_end_name,
    bool forced) {
    
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
        
        if(!forced) {
            get_aligned_motif_state(parent->data()->cur_state->end_states()[p],
                                    n_data->cur_state,
                                    n_data->ref_state);
        }
        //TODO need to comment this out of occasionally will not allow conversion from mt
        //if(sterics_ && _steric_clash(n_data)) { continue; }
    
        return tree_.add_data(n_data, (int)state->end_states().size(), parent->index(), p);
        
    }
    
    return -1;
}

int
MotifStateTree::add_mst(
    MotifStateTreeOP const & mst,
    int parent_index,
    int parent_end_index,
    String parent_end_name,
    bool forced) {
    
    int i = -1;
    int j = 0;
    auto index_dict = std::map<int, int>();
    
    for(auto const & n : *mst) {
        i++;
        
        
        if(i == 0) {
            j = add_state(n->data()->ref_state, parent_index, parent_end_index, parent_end_name,
                          forced);
        }
        else {
            int ind = index_dict[n->parent_index()];
            int pei = n->parent_end_index();
            j = add_state(n->data()->ref_state, ind, parent_end_index, "", forced);
            
        }
        
        index_dict[n->index()] = j;
        if(j == -1) {
            throw MotifStateTreeException("could not add motif state tree to this tree");
        }
    }
    
    return j;
    
    
}

int
MotifStateTree::setup_from_mt(
    MotifTreeOP const & mt) {
    
    int i = -1, j = -1;
    int parent_index = 1000, parent_end_index = -1;
    for(auto const & n : *mt) {
        i++;
        auto ms = RM::instance().motif_state(n->data()->name(),
                                             n->data()->end_ids()[0],
                                             n->data()->ends()[0]->name());
        if(i == 0) {
            add_state(ms);
        }
        else {
            parent_index = n->parent_index();
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

MotifTreeOP
MotifStateTree::to_motif_tree() {
    
    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", options_.get_bool("sterics"));
    int i = -1, j = -1;
    int parent_index = -1, parent_end_index = -1;
    for(auto const & n : tree_) {
        i++;
        MotifOP m;
        if(n->data()->ref_state->name() != "") {
            m = RM::instance().motif(n->data()->ref_state->name(), n->data()->ref_state->end_ids()[0]);
        }
        else {
            m = RM::instance().motif("", n->data()->ref_state->end_ids()[0]);
        }
        
        if(i == 0) {
            align_motif(n->data()->cur_state->end_states()[0],
                        m->ends()[0], m);
            mt->add_motif(m);
            continue;
        }
        
        parent_index = n->parent_index();
        parent_end_index = n->parent_end_index();
        if(parent_end_index == -1) {
            throw MotifStateTreeException("cannot convert to motif tree");
        }
        
        j = mt->add_motif(m, parent_index, parent_end_index);
        if(j == -1) {
            throw MotifStateTreeException("failed to add motif in to_motif_tree");
        }
        
        
    }
    
    return mt;
}

void
MotifStateTree::replace_state(
    int i,
    MotifStateOP const & new_state) {

    auto n = tree_.get_node(i);
    if(new_state->end_states().size() != n->data()->ref_state->end_states().size()) {
        throw MotifStateTreeException("attempted to replace a state with a different number of ends");
    }
    
    auto old_state = n->data()->ref_state;
    n->data()->ref_state = new_state;
    n->data()->cur_state = std::make_shared<MotifState>(new_state->copy());
    
    queue_.push(n);
    MotifStateTreeNodeOP current, parent;
    int pei;
    while(!queue_.empty()) {
        current = queue_.front();
        queue_.pop();
    
        parent = current->parent();
        if(parent == nullptr) { continue; }
        pei = current->parent_end_index();
        aligner_.get_aligned_motif_state(parent->data()->cur_state->end_states()[pei],
                                         current->data()->cur_state,
                                         current->data()->ref_state);
        
        for(auto const & c : current->children()) {
            if(c != nullptr) { queue_.push(c); }
        }
        
    }
}

String
MotifStateTree::topology_to_str() {
    String s;
    
    for(auto const & n : tree_) {
        s += n->data()->ref_state->name() + "," + n->data()->ref_state->end_names()[0] + ",";
        s += n->data()->ref_state->end_ids()[0] + "," + std::to_string(n->parent_index()) + ",";
        s += std::to_string(n->parent_end_index()) + " ";
    }
    s += "|";
    for(auto const & c : connections_) {
        s += c->to_str() + " ";
    }
    
    return s;
}

















