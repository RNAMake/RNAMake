//
//  motif_state_ensemble_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_ensemble_tree__
#define __RNAMake__motif_state_ensemble_tree__

#include <stdio.h>

//RNAMAke Headers

#include "data_structure/tree/tree.h"
#include "data_structure/tree/tree_node.h"
#include "motif/motif_state_ensemble.h"
#include "motif/motif_tree.h"
#include "motif_data_structures/motif_state_tree.h"

typedef TreeNodeOP<MotifStateEnsembleOP> MotifStateEnsembleTreeNodeOP;

class MotifStateEnsembleTree {
public:
    MotifStateEnsembleTree():
    tree_( TreeStatic<MotifStateEnsembleOP>()){}
    
    ~MotifStateEnsembleTree() {}
    
public:
    
    void
    setup_from_mt(
        MotifTreeOP const &);
    
    int
    add_ensemble(
        MotifStateEnsembleOP const & ensemble,
        int parent_index = -1,
        int parent_end_index = -1);
    
    MotifStateTreeOP
    to_mst();
    
    
public:

    size_t
    size() { return tree_.size(); }
    
    MotifStateEnsembleTreeNodeOP const &
    get_node(
        int i) {
        return tree_.get_node(i);
    }
    
public:
    inline
    MotifStateEnsembleTreeNodeOP const &
    last_node() { return tree_.last_node(); }
    

private:
    TreeStatic<MotifStateEnsembleOP> tree_;

};

typedef std::shared_ptr<MotifStateEnsembleTree> MotifStateEnsembleTreeOP;


#endif /* defined(__RNAMake__motif_state_ensemble_tree__) */
