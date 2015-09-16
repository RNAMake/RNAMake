//
//  motif_state_tree.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_tree__
#define __RNAMake__motif_state_tree__

#include <stdio.h>
#include <queue>

//RNAMake Headers
#include "base/types.h"
#include "base/base.h"
#include "data_structure/tree/tree.h"
#include "data_structure/tree/tree_node.h"
#include "motif/motif_state.h"
#include "motif/motif_state_aligner.h"
#include "motif/motif_tree.h"


class MotifStateTreeException : public std::runtime_error {
public:
    MotifStateTreeException(
        String const & message):
    std::runtime_error(message)
    {}
};

struct MSTNodeData {
    inline
    MSTNodeData(
        MotifStateOP const & nref_state):
    ref_state(nref_state),
    cur_state(std::make_shared<MotifState>(nref_state->copy()))
    {}
    
    inline
    BasepairStateOP const &
    get_end_state(
        String const & name) {
        return cur_state->get_end_state(name);
    }
    
    MotifStateOP ref_state, cur_state;
    
};

typedef std::shared_ptr<MSTNodeData> MSTNodeDataOP;
typedef TreeNodeOP<MSTNodeDataOP> MotifStateTreeNodeOP;

class MotifStateTree : public Base {
public:
    
    MotifStateTree() {
        tree_ = TreeStatic<MSTNodeDataOP>();
        aligner_ = MotifStateAligner();
        queue_ = std::queue<MotifStateTreeNodeOP>();
        setup_options(); update_var_options();
    }
    
    ~MotifStateTree() {}
    
public: //iterators
    
    typedef typename TreeStatic<MSTNodeDataOP>::iterator iterator;
    typedef typename TreeStatic<MSTNodeDataOP>::const_iterator const_iterator;
    
    iterator begin() { return tree_.begin(); }
    iterator end()   { return tree_.end(); }
    
    const_iterator begin() const { return tree_.begin(); }
    const_iterator end()   const { return tree_.end(); }
    
public:
    
    int
    add_state(
        MotifStateOP const & state,
        int parent_index=-1,
        int parent_end_index=-1,
        String parent_end_name="");

    int
    setup_from_mt(
        MotifTreeOP const &);
    
    MotifTreeOP
    to_motif_tree();
    
    void
    replace_state(
        int i, 
        MotifStateOP const &);
    
public: //motif tree wrappers
    
    void
    write_pdbs(
        String const & fname = "nodes") {
        to_motif_tree()->write_pdbs(fname);
    }
    
public:
    size_t
    size() { return tree_.size(); }
    
    MotifStateTreeNodeOP const &
    last_node() { return tree_.last_node(); }
    
    MotifStateTreeNodeOP const &
    get_node(
        int i) {
        return tree_.get_node(i);
    }
    
private:
    
    inline
    int
    _steric_clash(
        MSTNodeDataOP const & new_data) {
        
        float dist;
        for(auto const & n : tree_) {
            for(auto const & b1 : n->data()->cur_state->beads()) {
                for(auto const & b2 : new_data->cur_state->beads()) {
                    dist = b1.distance(b2);
                    if (dist < clash_radius_) { return 1; }
                }
            }
        }
        return 0;
    }
    
    
private:
    void
    setup_options();
    
    void
    update_var_options();
    
private:
    TreeStatic<MSTNodeDataOP> tree_;
    std::queue<MotifStateTreeNodeOP> queue_;
    MotifStateAligner aligner_;
    int sterics_;
    float clash_radius_;
    
};

typedef std::shared_ptr<MotifStateTree> MotifStateTreeOP;

#endif /* defined(__RNAMake__motif_state_tree__) */
