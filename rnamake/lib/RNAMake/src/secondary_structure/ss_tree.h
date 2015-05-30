//
//  ss_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_tree__
#define __RNAMake__ss_tree__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "data_structure/tree/tree.h"
#include "secondary_structure/ss_tree_node.h"
#include "secondary_structure/ss_tree_node.fwd.h"

class SS_Tree {
public:
    SS_Tree(
        String const &,
        String const &);

public:
    inline
    int
    size() { return (int)tree_.size(); }
    
    inline
    Node<SS_TreeNode*>*
    get_node(int pos) { return tree_.get_node(pos); }
    
    inline
    SS_TreeNode*
    get_data(int pos) { return tree_.get_data(pos); }
    
    inline
    std::vector<Node<SS_TreeNode*>*> const &
    nodes() { return tree_.nodes(); }
    
private:
    void
    _build_tree();
    
    std::vector<SS_TreeNode*>
    _build_tree_level(
        int,
        int);
    
    SS_TreeNode*
    _assign_new_node(
        int,
        int);
    
    int
    _get_brack_pair(
        int);
    
    int
    _get_dot_bounds(
        int,
        int);

private:
    String ss_, seq_;
    Tree<SS_TreeNode*> tree_;
    
};

#endif /* defined(__RNAMake__ss_tree__) */
