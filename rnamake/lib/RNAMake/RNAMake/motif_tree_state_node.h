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
    MotifTreeStateNode(
        MotifTreeStateOP const & mts,
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
        states_ = BasepairStateOPs(mts->end_states());
        int i = 0;
        for (auto const & s : mts_->end_states()) {
            if( s == NULL) { states_[i] = NULL; }
            else           {
                states_[i] = BasepairStateOP(new BasepairState(s->copy()));
            }
            i++;
        }
        children_ = MotifTreeStateNodeOPs(states_.size());
 
    }
    
    ~MotifTreeStateNode() {
     //   std::cout << "made it" << std::endl;
    }
    
    MotifTreeStateNode
    copy();
    
public:
    
    Ints const
    available_ends();
    
    void
    add_child(
        MotifTreeStateNodeOP const &,
        int const &);
    
    void
    remove_child(int const &);
    
    int
    parent_end_index();
    
    BasepairStateOP const &
    parent_end();
    
    MotifTreeStateNodeOPs
    children();
    
    void
    replace_mts(MotifTreeStateOP const &);
    
    
    inline
    BasepairStateOPs
    active_states() const {
        BasepairStateOPs states;
        for (auto const & s : states_ ) {
            if(s != NULL) { states.push_back(s); }
        }
        return states;
    }
    
    inline
    void
    disconnect() {
        parent_ = NULL;
        for(int i = 0; i < children_.size(); i++) {
            children_[i] = NULL;
        }
        for(auto & state : states_) {
            state = NULL;
        }
    }
    
public: //getters
    
    inline
    int const &
    level() const { return level_; }
    
    inline
    BasepairStateOPs const &
    states() const { return states_; }
    
    inline
    MotifTreeStateOP const &
    mts() const { return mts_; }
    
    inline
    Points const &
    beads() const { return beads_; }
    
    inline
    MotifTreeStateNodeOP const &
    parent() { return parent_; }
    
    inline
    int const &
    index() const { return index_; }

    inline
    int const &
    size() const { return size_; }
    
    inline
    float const &
    ss_score() const { return ss_score_; }
    
public: //setters
    
    inline
    void
    beads(Points const & nbeads) { beads_ = nbeads; }
    
    inline
    void
    size(int const & nsize) { size_ = nsize; }
    
    inline
    void
    ss_score(int const & nss_score) { ss_score_ = nss_score; }
    

private:
    MotifTreeStateOP mts_;
    MotifTreeStateNodeOP parent_;
    MotifTreeStateNodeOPs children_;
    int index_, lib_type_, level_, size_;
    float score_, ss_score_;
    Points beads_;
    BasepairStateOPs states_;
    
};

#endif /* defined(__RNAMake__motif_tree_state_node__) */
