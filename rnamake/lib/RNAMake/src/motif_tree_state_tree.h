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
    ~MotifTreeStateTree() {
        for(auto const & n : nodes_) { n->disconnect(); }
    }

public:
    MotifTreeStateNodeOP
    add_state(MotifTreeStateOP const & mts,
              MotifTreeStateNodeOP const & cparent,
              int const cparent_end = -1);

    MotifTree
    to_motiftree() const;
    
    PoseOP
    to_pose() const;
    
    int
    replace_state(MotifTreeStateNodeOP const &,
                  MotifTreeStateOP const &);
    void
    to_pdb(String fname="mtst.pdb",
           int include_head = 0) {
        MotifTree mt = to_motiftree();
        mt.to_pdb(fname, include_head);
    }
    
    void
    remove_node(
        MotifTreeStateNodeOP node = NULL);
    
public:
    inline
    MotifTreeStateNodeOPs const &
    nodes() const { return nodes_; }
    
    inline
    MotifTreeStateNodeOP const &
    last_node() { return last_node_; }

    
public:
    inline
    void
    sterics(int nsterics) { sterics_ = nsterics; }
    
    inline
    void
    clash_radius(float nclash_radius) { clash_radius_ = nclash_radius; }
    
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

typedef std::shared_ptr<MotifTreeStateTree> MotifTreeStateTreeOP;


#endif /* defined(__RNAMake__motif_tree_state_tree__) */
