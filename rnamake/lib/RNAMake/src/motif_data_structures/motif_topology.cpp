
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
                
                m = ResourceManager::getInstance().get_motif(name);
            }
            
            else {
                m  = ResourceManager::getInstance().get_motif(name, "",
                                                              current->data()->ends()[0]->name());
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
                m = ResourceManager::getInstance().get_motif(name);
            }
            else {
                m  = ResourceManager::getInstance().get_motif(current->data()->name(),
                                                              "", c_end_name);
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



















