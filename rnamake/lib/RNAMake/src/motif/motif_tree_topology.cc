//
//  motif_tree_topology.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include <queue>
#include "motif/motif_tree_topology.h"


MotifTreeTopology::MotifTreeTopology(SS_Tree const & ss_tree) {

    std::map<int, int> index_map_;
    index_map_[-1] = -1;
    SS_Nodes nodes;
    for(auto const & n : ss_tree) {
        int parent_index = n->parent_index();
        if(parent_index == -1) { continue; }
        
        auto parent = n->parent();
        if(parent->data()->type() == SS_NodeData::SS_Type::SS_BP &&
           n->data()->type()      == SS_NodeData::SS_Type::SS_BP) {
            nodes = SS_Nodes({parent, n});
            //String seq = combine_seq(nodes);
            SeqSS seq_ss = ss_tree.seq_from_nodes(nodes);
            for(int i = 0; i < seq_ss.ss.size(); i++) {
                std::cout << seq_ss.seq[i] << " " << seq_ss.ss[i] << std::endl;
            }
            
            exit(0);
            //auto mt = new MotifTopology(MotifTopology::MotifType::BASEPAIR, seq, "((+))");
            //exit(0);
        }
    }
}
