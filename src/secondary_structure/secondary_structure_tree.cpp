//
//  secondary_structure_tree.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/14/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//


#include <queue>

//RNAMake Headers
#include "secondary_structure/secondary_structure_tree.h"


namespace secondary_structure {
    
int
SecondaryStructureTree::add_motif(
    MotifOP const & m,
    int parent_index,
    int parent_end_index) {
    
    auto parent = tree_.last_node();
    
    if(parent_index != -1) {
        parent = tree_.get_node(parent_index);
    }
    
    auto m_copy = std::make_shared<Motif>(*m);
    if(parent == nullptr) {
        return tree_.add_data(m_copy, (int)m_copy->ends().size(), -1, -1);
    }
    
    auto avail_pos = tree_.get_available_pos(parent, parent_end_index);
    if(avail_pos.size() == 0) {
        throw std::runtime_error("cannot add motif to tree no available pos");
    }
    return tree_.add_data(m_copy, (int)m_copy->ends().size(), parent->index(), avail_pos[0]);
    
    return -1;
    
}
    
struct TreePreNode {
    MotifOP m;
    int parent_index, parent_end_index;
    
    inline
    TreePreNode() {}
    
    inline
    TreePreNode(
        MotifOP const & nm,
        int nparent_index,
        int nparent_end_index):
    m(nm),
    parent_index(nparent_index),
    parent_end_index(nparent_end_index)
    {}
};

SecondaryStructureTreeOP
tree_from_pose(PoseOP const & p) {
    auto sst = std::make_shared<SecondaryStructureTree>();
    auto start_m = MotifOP(nullptr);
    auto found = 0;
    for(auto const & m : p->motifs()) {
        found = 0;
        for(auto const & m1 : p->motifs()) {
            if(m == m1) { continue; }
            for(auto const & e : m1->ends()) {
                if(e == m->ends()[0]) {
                    found = 1;
                    break;
                }
            }
        }
        if(found == 0) {
            start_m = m;
            break;
        }
    }
    
    if(start_m == nullptr) {
        std::runtime_error("could not find start motif to build tree in secondary_structure::tree_from_pose");
    }
    
    auto open_nodes = std::queue<TreePreNode>();
    auto seen = MotifOPs();
    auto current = TreePreNode();
    open_nodes.push(TreePreNode{ start_m, -1 , -1});
    int i = 0;
    
    while(! open_nodes.empty()) {
        current = open_nodes.front();
        open_nodes.pop();
        
        seen.push_back(current.m);
        
        int pos = sst->add_motif(current.m, current.parent_index, current.parent_end_index);
        for(auto const & m : p->motifs()) {
            if(std::find(seen.begin(), seen.end(), m) != seen.end()) {
                continue;
            }
            i = -1;
            for(auto const & e : current.m->ends()){
                i++;
                if(i == 0) { continue; }
                if(e == m->ends()[0]) {
                    open_nodes.push(TreePreNode(m, pos, i));
                }
            }
        }
    }
    
    return sst;
}
    
    
}





















