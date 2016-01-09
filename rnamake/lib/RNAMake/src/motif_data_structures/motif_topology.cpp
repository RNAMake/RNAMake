//
//  motif_toplogy.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/8/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include <queue>

#include "motif_data_structures/motif_topology.h"
#include "resources/resource_manager.h"


MotifTreeOP
graph_to_tree(
    MotifGraphOP const & mg,
    GraphNodeOP<MotifOP> start,
    String last_end) {
    
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
    
    int highest = -1;
    auto parent = GraphNodeOP<MotifOP>(nullptr);
    auto p_end_name = String("");
    auto p_end_index = -1;
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
            continue;
        }
        
        highest = -1; p_end_index = -1; c_end_index = -1;
        p_end_name = ""; c_end_name = "";
        
        for(auto const & c : current->connections()) {
            if(c == nullptr) { continue; }
        }
        
    }
    
    return mt;
    
}



















