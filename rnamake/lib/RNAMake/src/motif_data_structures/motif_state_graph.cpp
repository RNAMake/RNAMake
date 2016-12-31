//
//  motif_state_graph.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/21/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "base/util.hpp"
#include "resources/resource_manager.h"
#include "motif_data_structures/motif_state_graph.hpp"


MotifStateGraph::MotifStateGraph():
    graph_(GraphStatic<MSNodeDataOP>()),
    clash_radius_(2.5),
    sterics_(1),
    options_(Options()),
    update_align_list_(1),
    align_list_(GraphNodeOPs<MSNodeDataOP>()),
    aligned_(std::map<int, int>()) { setup_options(); }

MotifStateGraph::MotifStateGraph(MotifGraphOP const & mg) : MotifStateGraph() {
    _setup_from_mg(mg);
}

MotifStateGraph::MotifStateGraph(MotifStateGraph const & msg) : MotifStateGraph() {
        
    graph_ = GraphStatic<MSNodeDataOP>(msg.graph_);
    // dear god this is horrible but cant figure out a better way to do a copy
    for(auto const & n : msg.graph_.nodes()) {
        graph_.get_node(n->index())->data() = std::make_shared<MSNodeData>(*n->data());
    }
}


void
MotifStateGraph::_setup_from_mg(MotifGraphOP const & mg) {
    aligned_ = mg->aligned();
    auto max_index = 0;
    for(auto const & n : *mg) {
        auto n_data = std::make_shared<MSNodeData>(n->data()->get_state());
        graph_.add_data(n_data, -1, -1, -1, n->data()->ends().size(), 1, n->index());
        if(n->index() > max_index) {
            max_index = n->index();
        }
    }
    graph_.index(max_index+1);
    for(auto const & c : mg->connections()) {
        graph_.connect(c->node_1()->index(), c->node_2()->index(),
                       c->end_index_1(), c->end_index_2());
    }
}


//add function helpers /////////////////////////////////////////////////////////////////////////////


GraphNodeOP<MSNodeDataOP>
MotifStateGraph::_get_parent(
    String const & m_name,
    int parent_index) {
    
    auto parent = graph_.last_node();
    
    //catch non existant parent
    try {
        if(parent_index != -1) { parent = graph_.get_node(parent_index); }
    }
    catch(GraphException const & e) {
        throw MotifStateGraphException(
            "could not add state: " + m_name + " with parent index: " +
            std::to_string(parent_index) + "there is no node with that index");
    }
    
    return parent;
}

Ints
MotifStateGraph::_get_available_parent_end_pos(
    GraphNodeOP<MSNodeDataOP> const & parent,
    int parent_end_index) {
    
    auto avail_pos = Ints();
    
    if(parent_end_index != -1) {
        int avail = parent->available_pos(parent_end_index);
        if(!avail) {
            throw MotifStateGraphException(
                    "cannot add state to graph as the end with index: " +
                    std::to_string(parent_end_index) +" you are trying to "
                    "add it to is already filled or does not exist");
        }
        
        if(parent_end_index == parent->data()->block_end_add()) {
            throw MotifStateGraphException(
                    "cannot add motif: to tree as the parent_end_index"
                    " supplied is blocked see class Motif");
        }
        
        avail_pos.push_back(parent_end_index);
    }
    
    else {
        auto avail_pos_temp = parent->available_children_pos();
        for(auto const & p : avail_pos_temp) {
            if(p == parent->data()->block_end_add()) { continue; }
            avail_pos.push_back(p);
        }
    }
    
    return avail_pos;
    
}

int
MotifStateGraph::_get_parent_index_from_name(
    GraphNodeOP<MSNodeDataOP> const & parent,
    String const & parent_end_name) {
    
    auto parent_end_index = -1;
    
    //called without selecting a parent and no motifs in graph
    if(parent == nullptr) {
        throw MotifStateGraphException(
            "cannot find parent_end_name in parent as none was specified");
    }
    
    
    try{
        parent_end_index = parent->data()->get_end_index(parent_end_name);
    }
    catch(MotifStateException) {
        throw MotifStateGraphException(
            "cannot find parent_end_name: " + parent_end_name + " while trying to "
            "add a state to graph");
    }
    
    int avail = parent->available_pos(parent_end_index);
    if(!avail) {
        throw MotifStateGraphException(
            "cannot add state to graph as the end with name: " +
            parent_end_name +" you are trying to "
            "add it to is already filled or does not exist");
    }
    
    if(parent_end_index == parent->data()->block_end_add()) {
        throw MotifStateGraphException(
            "cannot add state: to tree as the parent_name_index"
            " supplied is blocked see class Motif");
    }
    
    return parent_end_index;
}

int
MotifStateGraph::_get_connection_end(
    GraphNodeOP<MSNodeDataOP> const & node,
    String const & bp_name) {
    
    int node_end_index = -1;
    
    if(bp_name != "") {
        auto ei = node->data()->get_end_index(bp_name);
        
        if(!node->available_pos(ei)) {
            throw MotifStateGraphException(
                "cannot add connection with " + std::to_string(node->index()) + " and "
                "end name " + bp_name + " as the end is blocked");
        }
        
        node_end_index = ei;
        
    }
    
    else {
        auto node_indexes = node->available_children_pos();
        
        if(node_indexes.size() > 1) {
            throw MotifStateGraphException(
                "cannot connect nodes " + std::to_string(node->index()) + " its unclear "
                " which ends to attach as there is more then one possibility");
        }
        
        if(node_indexes.size() == 0) {
            throw MotifStateGraphException(
                "cannot connect nodes " + std::to_string(node->index())  + " there are "
                "no ends free ends to attach too");
        }
     
        node_end_index = node_indexes[0];
        
    }
    
    return node_end_index;
    
}


//add function  ////////////////////////////////////////////////////////////////////////////////////


int
MotifStateGraph::add_state(
    MotifStateOP const & state,
    int parent_index,
    int parent_end_index) {
    
    for(auto const & n : graph_.nodes()) {
        if(n->data()->uuid() == state->uuid()) {
            throw MotifStateGraphException(
                "cannot add state: " + state->name() + " to tree as its uuid is " +
                "already present in the tree");
        }
    }
    
    auto parent = _get_parent(state->name(), parent_index);
    
    if(parent == nullptr) {
        auto n_data = std::make_shared<MSNodeData>(state);
        auto pos =  graph_.add_data(n_data, -1, -1, -1, (int)state->end_states().size(), 1);
        if(pos != -1) {
            update_align_list_ = 1;
            aligned_[pos] = 0;
        }
        return pos;
    }
    
    auto avail_pos = _get_available_parent_end_pos(parent, parent_end_index);
    
    for(auto const & p : avail_pos) {
        auto n_data = std::make_shared<MSNodeData>(state);
        
        get_aligned_motif_state(parent->data()->cur_state->end_states()[p],
                                n_data->cur_state,
                                n_data->ref_state);
        
        if(sterics_ && _steric_clash(n_data)) { continue; }
        
        auto pos = graph_.add_data(n_data, parent->index(), p, 0, (int)state->end_states().size());
        if(pos != -1) {
            update_align_list_ = 1;
            aligned_[pos] = 1;
        }
        return pos;
    }
    
    return -1;
}

int
MotifStateGraph::add_state(
    MotifStateOP const & state,
    int parent_index,
    String const & parent_end_name) {
    
    auto parent = _get_parent(state->name(), parent_index);
    auto parent_end_index = _get_parent_index_from_name(parent, parent_end_name);
    return add_state(state, parent_index, parent_end_index);
}

void
MotifStateGraph::add_connection(
    int i,
    int j,
    String const & i_bp_name,
    String const & j_bp_name) {
    
    auto node_i = GraphNodeOP<MSNodeDataOP>(nullptr);
    auto node_j = GraphNodeOP<MSNodeDataOP>(nullptr);
    
    try {  node_i = graph_.get_node(i); }
    catch(GraphException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(i) +" does not exist");
    }
    
    try {  node_j = graph_.get_node(j); }
    catch(GraphException) {
        throw MotifTreeException(
            "cannot connect: " + std::to_string(i) + " " + std::to_string(j) + " as node " +
            std::to_string(j) +" does not exist");
    }
    
    auto node_i_ei = _get_connection_end(node_i, i_bp_name);
    auto node_j_ei = _get_connection_end(node_j, j_bp_name);
    
    auto node_i_end_name = node_i->data()->end_name(node_i_ei);
    auto node_j_end_name = node_j->data()->end_name(node_j_ei);
    
    graph_.connect(i, j, node_i_ei, node_j_ei);
    
}

void
MotifStateGraph::replace_state(
    int i,
    MotifStateOP const & new_state) {
    
    auto n = graph_.get_node(i);
    if(new_state->end_states().size() != n->data()->ref_state->end_states().size()) {
        throw MotifStateGraphException(
            "attempted to replace a state with a different number of ends");
    }
    
    auto old_state = n->data()->ref_state;
    n->data()->ref_state = new_state;
    n->data()->cur_state = std::make_shared<MotifState>(*new_state);
    n->data()->uuid(old_state->uuid());
    _align_states(i);
    
}


//remove functions /////////////////////////////////////////////////////////////////////////////////


void
MotifStateGraph::remove_state(int pos) {
    auto n = graph_.get_node(pos);
    graph_.remove_node(pos);
    aligned_.erase(pos);
    update_align_list_ = 1;
}

void
MotifStateGraph::remove_level(int level) {
    int pos = 0;
    while(pos < graph_.nodes().size()) {
        auto n = graph_.nodes()[pos];
        if(n->level() >= level) {
            remove_state(n->index());
            continue;
        }
        pos++;
    }
}


//motif graph wrappers /////////////////////////////////////////////////////////////////////////////

MotifGraphOP
MotifStateGraph::to_motif_graph() {
    auto mg = std::make_shared<MotifGraph>();
    mg->set_option_value("sterics", false);
    auto non_aligned_nodes = unaligned_nodes();
    auto seen_connections = std::map<GraphConnectionOP<MSNodeDataOP>, int>();
    auto index_hash = std::map<int, int>();
    auto j = 0;
    
    _update_align_list();
    for(auto const & n : align_list_) {
        auto m = RM::instance().motif(n->data()->name(), "", n->data()->end_name(0));
        
        if(std::find(non_aligned_nodes.begin(), non_aligned_nodes.end(), n) != non_aligned_nodes.end() ) {
            align_motif(n->data()->get_end_state(0), m->ends()[0], m);
            j = mg->add_motif(m, -1, -1, 1);
        }
        else {
            auto c = n->connections()[0];
            seen_connections[c] = 1;
            auto parent = c->partner(n->index());
            auto pei = c->end_index(parent->index());
            j = mg->add_motif(m, index_hash[parent->index()], pei);
        }
        
        index_hash[n->index()] = j;
        
    }
    
    for(auto const c : graph_.connections()) {
        if(seen_connections.find(c) != seen_connections.end()) { continue; }
        mg->add_connection(index_hash[c->node_1()->index()],
                           index_hash[c->node_2()->index()],
                           c->node_1()->data()->end_name(c->end_index_1()),
                           c->node_2()->data()->end_name(c->end_index_2()));
    }
    
    mg->update_indexes(index_hash);
    return mg;
}


//misc /////////////////////////////////////////////////////////////////////////////////////////////

void
MotifStateGraph::_update_align_list() {
    if(!update_align_list_) { return; }
    
    auto non_aligned_nodes = unaligned_nodes();
    auto open = std::queue<GraphNodeOP<MSNodeDataOP>>();
    auto used_nodes = std::map<GraphNodeOP<MSNodeDataOP>, int>();
    
    align_list_ = GraphNodeOPs<MSNodeDataOP>();
    
    for(auto const & start : non_aligned_nodes) {
        open.push(start);
        auto seen_nodes = std::map<GraphNodeOP<MSNodeDataOP>, int>();
        
        while (!open.empty()) {
            auto n = open.front();
            open.pop();
            
            seen_nodes[n] = 1;
            if(n->index() == start->index()) {
                align_list_.push_back(n);
            }
            else {
                // should this be block end?
                if(n->connections()[0] == nullptr) { continue; }
                auto c = n->connections()[0];
                auto parent = c->partner(n->index());
                if(used_nodes.find(parent) == used_nodes.end()) { continue; }
                align_list_.push_back(n);
            }
            
            used_nodes[n] = 1;
            int i = -1;
            for(auto const & c : n->connections()) {
                i++;
                if(i == n->data()->block_end_add() || c == nullptr) { continue; }
                auto partner_n = c->partner(n->index());
                if(seen_nodes.find(partner_n) != seen_nodes.end() ||
                   used_nodes.find(partner_n) != used_nodes.end()) { continue; }
                if(c->end_index(partner_n->index()) == partner_n->data()->block_end_add()) {
                    open.push(partner_n);
                }
                else if(partner_n->data()->cur_state->end_states().size() == 1) {
                    open.push(partner_n);
                }
            }
        }
    }
    
    update_align_list_ = 0;
}

void
MotifStateGraph::_align_states(int pos) {
    auto non_aligned_nodes = unaligned_nodes();
    
    int start = 1;
    if(pos != -1) { start = 0; }
    
    _update_align_list();
    for(auto const & n : align_list_) {
        if(start == 0) {
            if(n->index() == pos) { start = 1; }
            else { continue; }
        }
        
        if(element_in_vector(n, non_aligned_nodes)) { continue; }
        
        auto parent = n->connections()[0]->partner(n->index());
        auto pei = n->connections()[0]->end_index(parent->index());
        get_aligned_motif_state(parent->data()->get_end_state(pei),
                                n->data()->cur_state,
                                n->data()->ref_state);
    }
    
}


// getters /////////////////////////////////////////////////////////////////////////////////////////


GraphNodeOPs<MSNodeDataOP> const
MotifStateGraph::unaligned_nodes() const {
    auto nodes = GraphNodeOPs<MSNodeDataOP>();
    for(auto const & kv : aligned_) {
        if(kv.second == 0) { nodes.push_back(get_node(kv.first)); }
    }
    return nodes;
}



//option functions /////////////////////////////////////////////////////////////////////////////////


void
MotifStateGraph::setup_options() {
    options_.add_option("sterics", true, OptionType::BOOL);
    options_.add_option("clash_radius", 2.9f, OptionType::FLOAT);
    options_.lock_option_adding();
    update_var_options();
}

void
MotifStateGraph::update_var_options() {
    sterics_              = options_.get_bool("sterics");
    clash_radius_         = options_.get_int("clash_radius");
}



















