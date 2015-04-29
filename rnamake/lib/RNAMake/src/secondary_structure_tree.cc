//
//  secondary_structure_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "secondary_structure_tree.h"
#include <queue>

int
get_brack_pair(
    int pos,
    String const & ss) {
    int bracket_count = 0;
    int i = -1;
    for(auto const & s : ss) {
        i++;
        if( i < pos) { continue; }
        if (s == '(') { bracket_count++; }
        if (s == ')') { bracket_count--; }
        if(bracket_count == 0) { return i; }
    }
    throw "could not find bracket pair";
    
}

int
get_dot_bounds(
    int pos,
    String const & ss,
    int reverse) {
    
    if(!reverse) {
        for(int i = pos+1; i < ss.length(); i++) {
            if(ss[i] != '.') { return i-1; }
        }
        return (int)ss.length()-1;
    }
    else {
        for(int i = pos-1; i >= 0; i--) {
            if(ss[i] != '.') { return i+1; }
        }
        return 0;
    }
    return pos;
}


void
SecondaryStructureTree::_build_tree(
    String & ss,
    String & seq) {
    
    SecondaryStructureNodeOP node;
    if     ((ss[0] == '.' || ss.back() == '.') && nodes_.size() == 0) {
        node = SecondaryStructureNodeOP(new SSN_Bulge(ss, seq, 0, (int)ss.length()-1 , NULL));
    }
    
    else if(ss[0] == '(') {
        int pair = get_brack_pair(0, ss);
        node = SecondaryStructureNodeOP(new SSN_Basepair(ss, seq, 0, pair, NULL));
    }
    
    std::queue<SecondaryStructureNodeOP> open_nodes;
    SecondaryStructureNodeOP current;
    open_nodes.push(node);
    while (! open_nodes.empty()) {
        current = open_nodes.front();
        open_nodes.pop();
        current->assign_all_children(ss, seq);
        for(auto const & c : current->children()) { open_nodes.push(c); }
        nodes_.push_back(current);
    }
    
    return;
    
    int nway = 0;
    SecondaryStructureNodeOPs nodes_to_erase;
    for (auto & n : nodes_) {
        nway = 0;
        if(n->ss_type() != SSN_TWOWAY) { continue; }
        for(auto const & c : n->children()) {
            if(c->ss_type() != SSN_TWOWAY) { continue; }
            nodes_to_erase.push_back(c);
            nway=1;
        }
        if(nway) {
            SecondaryStructureNodeOP nway_node (new SSN_Junction(n));
            SecondaryStructureNodeOP parent = nway_node->parent();
            parent->replace_child(n, nway_node);
            nodes_to_erase.push_back(n);
            nodes_.push_back(nway_node);
        }
    }
    
    for(auto const & n : nodes_to_erase) {
        int pos = (int)(std::find(nodes_.begin(), nodes_.end(), n) - nodes_.begin());
        nodes_.erase(nodes_.begin()+pos);
    }
    
    //std::cout << nodes_.size() << std::endl;
    
}







