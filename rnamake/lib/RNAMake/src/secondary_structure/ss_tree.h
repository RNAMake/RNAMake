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
#include "secondary_structure/ss_tree_node.fwd.h"

class SS_Tree {
public:
    SS_Tree(
        String const &,
        String const &);

public:
    inline
    int
    size() { return (int)nodes_.size(); }
    
private:
    void
    _build_tree();
    
    SS_TreeNodeProtoOPs
    _build_tree_level(
        SS_TreeNodeProtoOP const &,
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
    SS_TreeNodeOPs nodes_;
    
};

#endif /* defined(__RNAMake__ss_tree__) */
