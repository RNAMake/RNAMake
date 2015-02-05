//
//  motif_tree_state_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/5/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_node__
#define __RNAMake__motif_tree_state_node__

#include <stdio.h>
#include "motif_tree_state.h"
#include "motif_tree_state_node.fwd.h"

class MotifTreeStateNode {
public:
    inline
    MotifTreeStateNode(
        MotifTreeState const & mts,
        int const & index,
        MotifTreeStateNodeOP const & parent,
        int const & lib_type):
        mts_(mts),
        index_(index),
        parent_(parent),
        lib_type_(lib_type)
    {
            
        if(parent == NULL) { level_ = 0; }
        else               { level_ = parent->level() + 1; }
        
        score_ = 1000;
        size_ = 0;
        ss_score_ = 0;
        beads_ = Points();
        states_ = BasepairStateOPs(mts.end_states());
        int i = 0;
        for (auto const & s : mts_.end_states()) {
            if( s == NULL) { states_[i] = NULL; }
            else           {
                states_[i] = BasepairStateOP(new BasepairState(s->copy()));
            }
            i++;
        }
        children_ = MotifTreeStateNodeOPs(states_.size());
 
    }
    
    ~MotifTreeStateNode() {}
    
    MotifTreeStateNode
    copy();
    
public:
    
    BasepairStateOPs const
    available_ends();
    
    void
    add_child(
        MotifTreeStateNodeOP const &,
        BasepairStateOP const &);
    
    int
    parent_end_index();
    
    BasepairStateOP const &
    parent_end();
    
    void
    replace_mts(MotifTreeState const &);
    
public: //getters
    
    inline
    int const &
    level() const { return level_; }
    
    inline
    BasepairStateOPs const &
    states() const { return states_; }
    
    inline
    MotifTreeState const &
    mts() const { return mts_; }

public: //setters
    
    inline
    void
    beads(Points const & nbeads) { beads_ = nbeads; }
    

private:
    MotifTreeState mts_;
    MotifTreeStateNodeOP parent_;
    MotifTreeStateNodeOPs children_;
    int index_, lib_type_, level_, size_;
    float score_, ss_score_;
    Points beads_;
    BasepairStateOPs states_;
    
};

#endif /* defined(__RNAMake__motif_tree_state_node__) */
