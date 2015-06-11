//
//  motif_tree_topology.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <queue>
#include <algorithm>

//RNAMake Headers
#include "motif/motif_tree_topology.h"

MotifTreeTopology::MotifTreeTopology(SS_Tree const & ss_tree) {
    
    bp_step_size_ = 2;
    graph_ = GraphStatic<MotifTopologyOP>();
    
    std::map<int, int> index_map_;
    int index = 0;
    
    for(auto const & n : ss_tree) {
        
        if(n->data()->type() == SS_NodeData::SS_Type::SS_BP) {
            auto nodes = get_bp_nodes(n, index_map_);
            if(nodes->size() < bp_step_size_) { continue; }
            auto mt_node = build_motif_topology_node(nodes, ss_tree, MotifTopologyType::BP_STEP);
            
            index += 1;
            for(auto const & node : *nodes) { index_map_[index] = node->index(); }

            if(graph_.size() == 0) {
                graph_.add_data(mt_node, -1, -1, -1, (int)mt_node->ends.size());
                continue;
            }
            
            auto ci = get_connectivity_info(mt_node);
            graph_.add_data(mt_node, ci.n_index, ci.i, ci.j, (int)mt_node->ends.size());
            
        }
        continue;
    }
}


ConnectivityInfo
MotifTreeTopology::get_connectivity_info(
    MotifTopologyOP const & mt_node) {

    ConnectivityInfo ci;
    
    int i = 0, j = 0;
    for(auto const & n : graph_) {
        i = 0;
        for(auto const & end_1 : n->data()->ends) {
            j = 0;
            for(auto const & end_2 : mt_node->ends) {
                if(end_1 == end_2) {
                    ci.i = i; ci.j = j; ci.n_index = n->index();
                    break;
                }
                j++;
            }
            i++;
        }
    }
    
    return ci;
    
}



VectorUP<SS_Node>
MotifTreeTopology::get_bp_nodes(
    SS_Node const & node,
    std::map<int, int> const & seen) {
    
    SS_Node current = node;
    auto nodes = std::make_unique<std::vector<SS_Node>>();
    while(current->data()->type() == SS_NodeData::SS_Type::SS_BP) {
        
        nodes->push_back(current);
        if(seen.find(current->index()) != seen.end()) { break; }
        
        current = current->parent();
        if(current == nullptr) { break; }
        
    }
    
    std::reverse(nodes->begin(), nodes->end());
    return nodes;
    
}

MotifTopologyOP
MotifTreeTopology::build_motif_topology_node(
    VectorUP<SS_Node> const & nodes,
    SS_Tree const & ss_tree,
    MotifTopologyType const & type) {
    
    auto sub_tree_data = ss_tree.seq_from_nodes(*nodes);
    
    std::vector<MotifTopologyBP> bps;
    int i = 0, j = 0;
    for(auto const & n : *nodes) {
        if(n->data()->type() != SS_NodeData::SS_Type::SS_BP) { continue; }
        i = 0;
        Strings chains;
        for(auto const & nb : n->data()->bounds()) {
            j = 0;
            for(auto const & b : sub_tree_data->bounds) {
                if(nb[0] == b[0]) {
                    chains.push_back(sub_tree_data->seq[j]);
                    break;
                }
                else if(nb[0] == b[1]) {
                    String temp = sub_tree_data->seq[j];
                    std::reverse(temp.begin(), temp.end());
                    chains.push_back(temp);
                    break;
                }
                
                j++;
                
            }
            i++;
        }
        
        if(chains.size() != 2) {
            throw std::runtime_error("did not find the correct number of chains for basepair");
        }
        
        auto bp = MotifTopologyBP(n->data()->seq(), chains, Ints({n->data()->bounds(0)[0], n->data()->bounds(1)[1]}));
        bps.push_back(bp);
    }
    
    auto mt_node = std::make_shared<MotifTopology>(type, sub_tree_data->seq, sub_tree_data->ss, bps);
    return mt_node;
    
}








