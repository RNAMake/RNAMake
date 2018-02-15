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
#include "motif_data_structures/motif_tree.h"
#include "motif_data_structures/motif_state_tree.h"
#include "motif_data_structures/motif_connection.h"
#include "math/euler.h"

class MotifStateTreeEnsembleException : public std::runtime_error {
public:
    MotifStateTreeEnsembleException(
        String const & message):
    std::runtime_error(message)
    {}
};

typedef TreeNodeOP<MotifStateEnsembleOP> MotifStateEnsembleTreeNodeOP;

/*const double _EPS = 2.22044604925e-16 * 4.0;

//assumes 3x3 matrices
inline
void
calc_euler2(
    Matrix & M,
    Vector & euler) {
    
    double cy = sqrt(M.xx()*M.xx() + M.yx()*M.yx());
    if(cy > _EPS) {
        euler[0] = atan2( M.zy(), M.zz());
        euler[1] = atan2(-M.zx(), cy);
        euler[2] = atan2( M.yx(), M.xx());
    }
    else {
        euler[0] = atan2( M.yz(), M.yy());
        euler[1] = atan2(-M.zx(), cy);
        euler[2] = 0.0;
    }
    for(int i = 0; i < 3; i++){
        if(euler[i] > 6.14) {
            euler[i] -= 6.14;
        }
        if(euler[i] < 0) {
            euler[i] += 6.14;
        }
    }
    
    //'sxyz': (0, 0, 0, 0)
    //_NEXT_AXIS = [1, 2, 0, 1]
    
}*/


class MotifStateEnsembleTree {
public:
    
    MotifStateEnsembleTree();
    
    MotifStateEnsembleTree(
        MotifTreeOP const &);
    
    MotifStateEnsembleTree(
        MotifStateTreeOP const &);
    
    ~MotifStateEnsembleTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MotifStateEnsembleOP>::iterator iterator;
    typedef typename TreeStatic<MotifStateEnsembleOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public: // add functions
        
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
    
    inline
    MotifStateEnsembleTreeNodeOP const &
    last_node() { return tree_.last_node(); }
    

private:
    TreeStatic<MotifStateEnsembleOP> tree_;
    MotifConnections connections_;

};

typedef std::shared_ptr<MotifStateEnsembleTree> MotifStateEnsembleTreeOP;

class MotifStateEnsembleTreeEnumerator {
public:
    MotifStateEnsembleTreeEnumerator(
        MotifStateEnsembleTreeOP const & mtst):
    mtst_(mtst)
    {}
    
    ~MotifStateEnsembleTreeEnumerator() {}
    
public:
    
    void
    record(
        String fname="test");
    
public:
    MotifStateEnsembleTreeOP mtst_;
    
};


#endif /* defined(__RNAMake__motif_state_ensemble_tree__) */
