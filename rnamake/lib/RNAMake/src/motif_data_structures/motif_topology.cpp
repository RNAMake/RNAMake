
//  motif_toplogy.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <queue>
#include <map>

#include "motif_data_structures/motif_topology.h"
#include "resources/resource_manager.h"


MotifTreeOP
GraphtoTree::convert(
    MotifGraphOP const & mg,
    GraphNodeOP<MotifOP> start,
    int start_end_index,
    GraphNodeOP<MotifOP> last_node) {
    
    mt_ = std::make_shared<MotifTree>();
    mt_->set_option_value("sterics", false);
    
    auto start_node = _get_start_node(mg, start, start_end_index);
    auto last_node_to_add = _GraphtoTreeNodeOP(nullptr);
    auto current = _GraphtoTreeNodeOP(nullptr);
    auto open_nodes = std::queue<_GraphtoTreeNodeOP>();
    auto new_nodes = _GraphtoTreeNodeOPs();
    auto seen_nodes = std::map<int, int>();
    open_nodes.push(start_node);
    seen_nodes[start_node->node->index()] = 1;

    while( open_nodes.size() > 0) {
        current = open_nodes.front();
        open_nodes.pop();
        
    
        if(last_node == current->node) {
            last_node_to_add = current;
            continue;
        }
        
        if(current->parent == nullptr) {
            mt_->add_motif(current->motif);
        }
        else {
            auto parent = mt_->get_node(current->parent->data()->id());
            auto new_parent_index = parent->index();
            
            //currently a hack to make sure that if we start building not at an end and trying to add
            //somethign on the 0th end of a HELIX this stops it
            if(parent->data()->mtype() == MotifType::HELIX && parent->children()[1] != nullptr) {
                seen_nodes.erase(current->node->index());
                continue;
            }
            
            //can happen sometimes if you build from the middle of the graph
            if(current->parent_end_index == 0) {
                continue;
            }
            
            auto pos = mt_->add_motif(current->motif, new_parent_index, current->parent_end_index);
            if(pos == -1) {
                continue;
            }
        }
        
        new_nodes = _get_new_nodes(current);
        for(auto const & n : new_nodes) {
            if(n == last_node_to_add) { continue; }
            
            if(seen_nodes.find(n->node->index()) == seen_nodes.end()) {
                seen_nodes[n->node->index()] = 1;
                open_nodes.push(n);
            }
        }
        
    }
    
    if(last_node_to_add != nullptr) {
        auto new_parent_index = mt_->get_node(last_node_to_add->parent->data()->id())->index();
        mt_->add_motif(last_node_to_add->motif,
                       new_parent_index,
                       last_node_to_add->parent_end_index);
    }
    
    
    mt_->set_option_value("sterics", false);
    return mt_;
    
}

GraphtoTree::_GraphtoTreeNodeOP
GraphtoTree::_get_start_node(
    MotifGraphOP const & mg,
    GraphNodeOP<MotifOP> const & start,
    int start_end_index) {
    
    auto start_n = _GraphtoTreeNodeOP(nullptr);
    if(start == nullptr) {
        auto not_aligned = mg->unaligned_nodes();
        if(not_aligned.size() == 0) {
            throw MotifTopologyException("cannot convert graph to tree no starting point");
        }
        auto start = not_aligned[0];
        start_n = std::make_shared<_GraphtoTreeNode>(nullptr, 0, start, start->data());
    }
    else {
        start_n = std::make_shared<_GraphtoTreeNode>(nullptr, 0, start, start->data());
        start_n->motif = _get_reoriented_motif(start->data(), start_end_index);

    }
    
    return start_n;
    
}

MotifOP
GraphtoTree::_get_reoriented_motif(
    MotifOP const & m,
    int end_index) {
    
    // dont need to do anything already oriented
    if(end_index == 0) { return m; }
    
    if(m->mtype() != MotifType::HELIX) {
        try {
            auto new_m = RM::instance().motif(m->name(), "", m->end_name(end_index));
            new_m->copy_uuids_from_motif(*m);
            return new_m;
        }
        catch(ResourceManagerException const & e) {
            throw MotifTopologyException(
                String("cannot convert graph to tree because ") + e.what());
        }
    }
    
    else {
        if(m->name().substr(0, 5) == "HELIX") {
            auto new_m = RM::instance().motif(m->name());
            new_m->copy_uuids_from_motif(*m);
            return new_m;

            
        }
        
        else {
            auto spl = split_str_by_delimiter(m->name(), "=");
            auto new_name = spl[1] + "=" + spl[0];
            auto new_m = RM::instance().motif(m->name());
            //since each basepair step has a different motif to represent it
            //these are likly two different motifs cannot do copy_uuids_from_motif
            new_m->id(m->id());
            return new_m;
        }
        
    }
    
}


GraphtoTree::_GraphtoTreeNodeOPs
GraphtoTree::_get_new_nodes(
    GraphtoTree::_GraphtoTreeNodeOP const & current) {
    
    auto new_nodes = _GraphtoTreeNodeOPs();
    auto partner = GraphNodeOP<MotifOP>();
    auto n = TreeNodeOP<MotifOP>();
    auto found = 1;
    
    for(auto const & c : current->node->connections()) {
        if(c == nullptr) { continue; }
        
        found = 1;
        partner = c->partner(current->node->index());
        // check to see if motif is already in tree
        try {
            n = mt_->get_node(partner->data()->id());
        } catch(...) { found = 0; }
        
        
        if(found) { continue; }
        
        auto parent_end_index = c->end_index(current->node->index());
        auto node_end_index = c->end_index(partner->index());
        auto new_n = std::make_shared<_GraphtoTreeNode>(current->node, parent_end_index,
                                                        partner, partner->data());
        new_n->motif = _get_reoriented_motif(new_n->motif, node_end_index);
        new_n->parent_end_index = _get_new_parent_end_index(current->node, c);
    
        //if(new_n->parent_end_index != 0) { continue; }
        
        new_nodes.push_back(new_n);

     }
    
    return new_nodes;
    
}

int
GraphtoTree::_get_new_parent_end_index(
    GraphNodeOP<MotifOP> const & parent,
    GraphConnectionOP<MotifOP> const & c) {
 
    if(parent->data()->mtype() != MotifType::HELIX) {
        auto parent_end_index = c->end_index(parent->index());
        auto parent_end_name = parent->data()->end_name(parent_end_index);
        auto tree_parent = mt_->get_node(parent->data()->id());
        
        int i = -1;
        for(auto const & end : tree_parent->data()->ends()) {
            i++;
            if(end->name() == parent_end_name) {
                return i;
            }
        }
        
        throw std::runtime_error("did not find original end something went really wrong");
    }
    
    //helices always go end 0 to 1
    else {
        return 1;
    }
    
}


MotifTreeOP
graph_to_tree(
    MotifGraphOP const & mg,
    GraphNodeOP<MotifOP> start,
    BasepairOP last_end) {
    
    if(start == nullptr) {
        start = mg->oldest_node();
    }
    
    auto mt = std::make_shared<MotifTree>();
    mt->set_option_value("sterics", false);
    auto open = std::queue<GraphNodeOP<MotifOP>>();
    auto current = GraphNodeOP<MotifOP>(nullptr);
    auto index = 0;
    auto m = MotifOP(nullptr);
    open.push(start);
    
    auto seen_nodes = std::map<GraphNodeOP<MotifOP>, int>();
    auto current_nodes = std::map<GraphNodeOP<MotifOP>, int>();
    auto seen_connections = std::map<String, int>();
    current_nodes[start] = 1;
    
    int highest = -1, pos = -1;
    auto p = GraphNodeOP<MotifOP>(nullptr);
    auto p_new = GraphNodeOP<MotifOP>(nullptr);
    auto p_end_name = String("");
    auto p_end_index = -1;
    auto p_index = -1;
    auto c_end_index = -1;
    auto c_end_name = String("");
    
    while( open.size() > 0) {
        current = open.front();
        open.pop();
        
        if(index == 0) {
            int free_end = -1;
            int i = 0;
            for(auto const & c : current->connections()) {
                if(c == nullptr) { free_end = i; break; }
            }
            assert(free_end != -1 && "no free end, cannot convert to motif tree");
            
            auto name = current->data()->name();
            if(name[2] == '=') {
                if(free_end == 1) {
                    auto spl = split_str_by_delimiter(name, "=");
                    name = spl[1] + "=" + spl[0];
                }
                
                m = RM::instance().motif(name);
            }
            
            else {
                m = RM::instance().motif(name, "", current->data()->ends()[0]->name());
            }
            mt->add_motif(m);
            mt->get_node(0)->data()->id(current->data()->id());
        }
        
        else {
        
            highest = -1; p_end_index = -1; c_end_index = -1;
            p_end_name = ""; c_end_name = "";
            p = nullptr;
            
            for(auto const & c : current->connections()) {
                if(c == nullptr) { continue; }
                p_new = c->partner(current->index());
                if(seen_nodes.find(p_new) == seen_nodes.end()) { continue; }
                p_index = seen_nodes[p_new];
                if(highest > p_index) { continue; }
                highest = p_index;
                p = p_new;
                p_end_index = c->end_index(p->index());
                p_end_name = p->data()->ends()[p_end_index]->name();
                c_end_index = c->end_index(current->index());
                c_end_name = current->data()->ends()[c_end_index]->name();
                
            }
            
            if(last_end != nullptr) {
                if(last_end->name() == p_end_name &&
                   seen_nodes.size() != mg->size()-1) {
                    open.push(current);
                    continue;
                }
            }
            
            
            assert(p != nullptr && "could not find parent something really wrong");
            //std::cout << p->data()->name() << " " << p_end_index << std::endl;
            if(p->data()->name()[2] == '=') {
                p_end_index = 1;
                auto act_parent = mt->get_node(seen_nodes[p]);
                p_end_name = act_parent->data()->ends()[1]->name();
            }
            
            if(current->data()->name()[2] == '=' && c_end_index == 1) {
                auto spl = split_str_by_delimiter(current->data()->name(), "=");
                auto name = spl[1] + "=" + spl[0];
                m = RM::instance().motif(name);
            }
            else {
                m = RM::instance().motif(current->data()->name(), "", c_end_name);
            }
            
            seen_connections[std::to_string(p_index) + " " + std::to_string(index)] = 1;
            pos = mt->add_motif(m, p_index, p_end_name);
            mt->get_node(pos)->data()->id(current->data()->id());
            assert(pos != -1 && "did not sucessfully add motif to tree during conversion");
            
        }
        
        seen_nodes[current] = index;
        current_nodes.erase(current);
        index++;
        
        for(auto const & c : current->connections()) {
            if(c == nullptr) { continue; }
            p = c->partner(current->index());
            if(seen_nodes.find(p) != seen_nodes.end()) { continue; }
            if(current_nodes.find(p) != current_nodes.end()) { continue; }
            open.push(p);
            current_nodes[p] = 1;
        }
    }
    
    String key1, key2;
    for(auto const & n : *mg) {
        for(auto const & c : n->connections()) {
            if(c == nullptr) { continue; }
            key1 = std::to_string(seen_nodes[c->node_1()]) + " " + std::to_string(seen_nodes[c->node_2()]);
            key2 = std::to_string(seen_nodes[c->node_2()]) + " " + std::to_string(seen_nodes[c->node_1()]);
            if(seen_connections.find(key1) == seen_connections.end() &&
               seen_connections.find(key2) == seen_connections.end()) {
                seen_connections[key1] = 1;
                mt->add_connection(seen_nodes[c->node_1()], seen_nodes[c->node_2()],
                                   "", "");
            }
        }
    }

    
    return mt;
    
}



















