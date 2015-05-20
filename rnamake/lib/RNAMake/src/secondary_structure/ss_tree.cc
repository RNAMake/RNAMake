//
//  ss_tree.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <queue>

#include "secondary_structure/ss_tree.h"
#include "secondary_structure/ss_tree_node.h"

SS_Tree::SS_Tree(
    String const & ss,
    String const & seq):
    ss_(ss),
    seq_(seq),
    nodes_(SS_TreeNodeOPs()) {
        
    _build_tree();
}

void
SS_Tree::_build_tree() {
    int xb = 0, yb = (int)seq_.size();
    
    SS_TreeNodeProtoOP current(nullptr, xb, yb);
    std::queue<SS_TreeNodeProtoOP> open_nodes;
    open_nodes.push(current);

    while(! open_nodes.empty()) {
        
        SS_TreeNodeProtoOPs next_level = _build_tree_level(current )
        
        
        if(ss_[xb] == '(') {
            seqs[0] = ss_[xb];
            seqs[1] = ss_[yb];
            current = SS_TreeNodeOP(new SS_TreeNodeBasepair(seqs, xb, yb));
            nodes_.push_back(current);
            xb += 1; yb -= 1;
            continue;
        }
    }
}

SS_TreeNodeProtoOPs
SS_Tree::_build_tree_level(
    SS_TreeNodeProtoOP const & current,
    int xb,
    int yb) {
    
    //while(xb, )
    
    
}


int
SS_Tree::_get_brack_pair(
    int pos) {
    
    int bracket_count = 0;
    int i = -1;
    for(auto const & s : ss_) {
        i++;
        if( i < pos) { continue; }
        if (s == '(') { bracket_count++; }
        if (s == ')') { bracket_count--; }
        if(bracket_count == 0) { return i; }
    }
    
    throw std::runtime_error("could not find bracket pair");

}

int
SS_Tree::_get_dot_bounds(
    int pos,
    int reverse) {

    if(!reverse) {
        for(int i = pos+1; i < ss_.length(); i++) {
            if(ss_[i] != '.') { return i-1; }
        }
        return (int)ss_.length()-1;
    }
    else {
        for(int i = pos-1; i >= 0; i--) {
            if(ss_[i] != '.') { return i+1; }
        }
        return 0;
    }
    return pos;
}


