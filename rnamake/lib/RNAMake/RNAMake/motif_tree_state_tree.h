//
//  motif_tree_state_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_tree__
#define __RNAMake__motif_tree_state_tree__

#include <stdio.h>
#include "motif_tree.h"
#include "motif_tree_state.h"
#include "motif_tree_state_node.h"
#include "motif_tree_state_node_aligner.h"


class MotifTreeStateTree {
public:
    MotifTreeStateTree();
    MotifTreeStateTree(MotifTreeStateOP const &);
    ~MotifTreeStateTree() {}

public:
    MotifTreeStateNodeOP
    add_state(MotifTreeStateOP const & mts,
              MotifTreeStateNodeOP const & cparent,
              int const cparent_end = -1);

    MotifTree
    to_motiftree();
    
    int
    replace_state(MotifTreeStateNodeOP const &,
                  MotifTreeStateOP const &);
    
public:
    inline
    MotifTreeStateNodeOPs const &
    nodes() { return nodes_; }
    
public:
    inline
    void
    sterics(int nsterics) { sterics_ = nsterics; }
    
private:
    int
    _steric_clash(
        MotifTreeStateNodeOP const &);

private:
    MotifTreeStateNodeOPs nodes_;
    MotifTreeStateNodeOP last_node_;
    MotifTreeStateNodeAligner aligner_;
    float clash_radius_;
    int sterics_;
    
    
};

MotifTreeState
ref_mts();


#endif /* defined(__RNAMake__motif_tree_state_tree__) */
