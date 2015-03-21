//
//  motif_tree_state_search_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/23/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_tree_state_search_node__
#define __RNAMake__motif_tree_state_search_node__

#include <stdio.h>
#include "types.h"
#include "motif_tree_state.h"
#include "motif_tree_state_node.h"
#include "motif_tree_state_search_node.fwd.h"

class MotifTreeStateSearchNode  {
public:
    MotifTreeStateSearchNode(
        MotifTreeStateOP const & ,
        MotifTreeStateSearchNodeOP const &,
        int);
    
    MotifTreeStateSearchNode (
        MotifTreeStateSearchNodeOP const & n):
    mts_ ( n->mts_ ),
    score_ ( n->score_ ),
    size_ ( n->size_ ),
    level_ ( n->level_ ),
    ss_score_ ( n->ss_score_ ),
    parent_ ( n->parent_ ),
    active_ ( n->active_ ),
    node_counts_ ( n->node_counts_ ),
    beads_ ( n->beads_ ),
    lib_type_ ( n->lib_type_ )
    
    {
        states_ = BasepairStateOPs(mts_->end_states().size());
        int i = 0;
        for(auto const & s : n->states()) {
            states_[i] = BasepairStateOP( new BasepairState(s->copy()));
            i++;
        }
    }
        
    ~MotifTreeStateSearchNode() {}
    
public:
    
    //MotifTreeStateSearchNodeOP
    //copy();
    
public:
    
    inline
    float ss_score() const { return ss_score_; }
    
    inline
    int size() const { return size_; }
    
    inline
    int level() const { return level_; }
    
    inline
    int lib_type() const { return lib_type_; }

    inline
    MotifTreeStateOP const &
    mts() const { return mts_; }
    
    inline
    MotifTreeStateSearchNodeOP const &
    parent() const { return parent_; }
    
    inline
    float const &
    score() const { return score_; }
    
    inline
    BasepairStateOPs const &
    states() const { return states_; }
    
    inline
    Points const &
    beads() const { return beads_; }

public:
    
    inline
    void
    beads(Points const & nbeads) { beads_ = nbeads; }
    
    inline
    void
    score(float const & nscore) { score_ = nscore; }
    
    inline
    void
    parent(MotifTreeStateSearchNodeOP const & nparent) { parent_ = nparent; }
    
    inline
    void
    level(int const & nlevel) { level_ = nlevel; }
    
    inline
    void
    lib_type(int const & nlib_type) { lib_type_ = nlib_type; }
    
    inline
    void
    ss_score(float const & nss_score) { ss_score_ = nss_score; }
    
    inline
    void
    size(int const nsize ) { size_ = nsize; }
    
public:
    
    inline
    int lib_type_usage(int lib_type) { return node_counts_[lib_type]; }
    
    inline
    Ints const &
    active() { return active_; }
    
    inline
    Ints const &
    node_counts() { return node_counts_; }
    
    void
    setup_node_count(int size) {
        node_counts_.resize(size);
        if(lib_type_ < 0) { return; }
        node_counts_[lib_type_]++;
    }
    
    void
    update_node_count() {
        node_counts_ = parent_->node_counts();
        if(lib_type_ < 0) { return; }
        node_counts_[lib_type_]++;
    }
    
    BasepairStateOPs
    active_states() const;
    
    void
    replace_mts(MotifTreeStateOP const & nmts);
    
    inline
    int
    steric_clash(int count = 9999) {
        int c = 0;
        float dist;
        MotifTreeStateSearchNodeOP current = parent_;
        while (current != NULL) {
            for (auto const b1 : beads_) {
                for (auto const b2 : current->beads()) {
                    dist = b1.distance(b2);
                    if (dist < 2.9) { return 1; }
                }
            }
            current = current->parent();
            c++;
            if (c > count) { break; }
        }
        return 0;
    }
    
    inline
    void
    update_stats() {
        size_ = parent_->size_ + mts_->size();
        ss_score_ = parent_->ss_score_ + mts_->score();
    }
    
private:
    MotifTreeStateOP mts_;
    MotifTreeStateSearchNodeOP parent_;
    Points beads_;
    BasepairStateOPs states_;
    Ints node_counts_, active_;
    int lib_type_, level_, size_;
    float ss_score_, score_;
    
};


#endif /* defined(__RNAMake__motif_tree_state_search_node__) */
