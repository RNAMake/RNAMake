//
//  motif_tree_topology.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>
#include "motif/motif_tree_topology.h"

/*
MotifTreeTopology::MotifTreeTopology(SS_Tree const & ss_tree) {

    std::map<int, int> index_map_;
    index_map_[-1] = -1;
    SS_Nodes nodes;
    for(auto const & n : ss_tree.nodes()) {
        int parent_index = n->parent_index();
        if(parent_index == -1) { continue; }
        
        auto parent = n->parent();
        if(parent->data()->type() == SS_TreeNode::SS_Type::SS_BP &&
           n->data()->type()      == SS_TreeNode::SS_Type::SS_BP) {
            nodes = SS_Nodes({parent, n});
            String seq = combine_seq(nodes);
            auto mt = new MotifTopology(MotifTopology::MotifType::BASEPAIR, seq, "((+))");
            exit(0);
        }
    }
}

String
MotifTreeTopology::combine_seq(
    SS_Nodes const & nodes) {
    
    String s1, s2;
    
    for(auto const & n : nodes) {
        Strings seqs = n->data()->seqs();
        s1 += seqs[0];
        s2 += seqs[1];
        if(seqs.size() != 2) {
            throw std::runtime_error("not implmented: cannot combine nodes with 3 strands");
        }
    }
    
    return s1 + "+" + s2;
    
}
*/