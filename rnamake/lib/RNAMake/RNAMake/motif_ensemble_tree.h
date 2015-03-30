//
//  motif_ensemble_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_ensemble_tree__
#define __RNAMake__motif_ensemble_tree__

#include <stdio.h>
#include "motif_ensemble.h"
#include "motif_ensemble_tree.fwd.h"
#include "motif_tree_state_tree.h"

class MotifEnsembleTree {
public:
    MotifEnsembleTree();
    MotifEnsembleTree(MotifEnsemble const &);
    MotifEnsembleTree(MotifTreeStateTree const &);
    ~MotifEnsembleTree() {}
    
public:
    MotifEnsembleTreeNodeOP
    add_ensemble(
        MotifEnsemble const & ensemble,
        MotifEnsembleTreeNodeOP parent = NULL,
        int parent_end_index = -1);
    
    MotifTreeStateTree
    get_mtst();
    
public:
    
    inline
    MotifEnsembleTreeNodeOPs const &
    nodes() { return nodes_; }
    

private:
    MotifEnsembleTreeNodeOPs nodes_;
    MotifEnsembleTreeNodeOP last_node_;
};

class MotifEnsembleTreeNode {
public:
    MotifEnsembleTreeNode(
        MotifEnsemble const & motif_ensemble,
        MotifEnsembleTreeNodeOP const & parent,
        int const index):
        motif_ensemble_(motif_ensemble),
        parent_(parent),
        index_(index)
    {
        children_ = MotifEnsembleTreeNodeOPs(motif_ensemble_.motif_states()[0].mts->end_states().size());
        start_index_ = motif_ensemble_.motif_states()[0].mts->start_index();
    }
    
    ~MotifEnsembleTreeNode() {}
    
public:
    
    Ints
    available_ends();
    
    void
    add_child(MotifEnsembleTreeNodeOP const &,
              int);
    
public: // getters
    inline
    MotifEnsemble const &
    motif_ensemble() const { return motif_ensemble_; }
    
    inline
    int const &
    index() const { return index_; }
    
    inline
    MotifEnsembleTreeNodeOPs const &
    children() { return children_; }
    
    inline
    MotifEnsembleTreeNodeOP const &
    parent() { return parent_; }

    inline
    int const
    parent_index() {
        MotifEnsembleTreeNodeOPs children = parent_->children();
        int i = -1;
        for(auto const & n : children) {
            i++;
            if(n == NULL) { continue; }
            if(n->index() == index_) { return i; }
        }
        return -1;
    }
    
private:
    MotifEnsemble motif_ensemble_;
    MotifEnsembleTreeNodeOP parent_;
    MotifEnsembleTreeNodeOPs children_;
    int index_;
    int start_index_;
};


#endif /* defined(__RNAMake__motif_ensemble_tree__) */
