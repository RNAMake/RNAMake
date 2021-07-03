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
#include "motif_data_structure/motif_tree.h"
#include "motif_data_structure/motif_state_tree.h"
#include "motif_data_structure/motif_connection.h"
#include "math/euler.h"

namespace motif_data_structure {

class MotifStateTreeEnsembleException : public std::runtime_error {
public:
    MotifStateTreeEnsembleException(
            String const & message) :
            std::runtime_error(message) {}
};

typedef data_structure::tree::TreeNodeOP<motif::MotifStateEnsembleOP> MotifStateEnsembleTreeNodeOP;

/*const double _EPS = 2.22044604925e-16 * 4.0;

//assumes 3x3 matrices
inline
void
math::calc_euler2(
    math::Matrix & M,
    math::Vector & euler) {
    
    double cy = sqrt(M.get_xx()*M.get_xx() + M.get_yx()*M.get_yx());
    if(cy > _EPS) {
        euler[0] = atan2( M.get_zy(), M.get_zz());
        euler[1] = atan2(-M.get_zx(), cy);
        euler[2] = atan2( M.get_yx(), M.get_xx());
    }
    else {
        euler[0] = atan2( M.get_yz(), M.get_yy());
        euler[1] = atan2(-M.get_zx(), cy);
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

    typedef typename data_structure::tree::TreeStatic<motif::MotifStateEnsembleOP>::iterator iterator;
    typedef typename data_structure::tree::TreeStatic<motif::MotifStateEnsembleOP>::const_iterator const_iterator;

    iterator begin() { return tree_.begin(); }

    iterator end() { return tree_.end(); }

    const_iterator begin() const { return tree_.begin(); }

    const_iterator end() const { return tree_.end(); }

public: // add functions

    int
    add_ensemble(
            motif::MotifStateEnsembleOP const & ensemble,
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
    data_structure::tree::TreeStatic<motif::MotifStateEnsembleOP> tree_;
    MotifConnections connections_;

};

typedef std::shared_ptr<MotifStateEnsembleTree> MotifStateEnsembleTreeOP;

class MotifStateEnsembleTreeEnumerator {
public:
    MotifStateEnsembleTreeEnumerator(
            MotifStateEnsembleTreeOP const & mtst) :
            mtst_(mtst) {}

    ~MotifStateEnsembleTreeEnumerator() {}

public:

    void
    record(
            String fname = "test");

public:
    MotifStateEnsembleTreeOP mtst_;

};

}


#endif /* defined(__RNAMake__motif_state_ensemble_tree__) */
