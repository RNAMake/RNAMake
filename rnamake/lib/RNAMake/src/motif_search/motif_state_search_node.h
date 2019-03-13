//
//  motif_state_search_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_state_search_node__
#define __RNAMake__motif_state_search_node__

#include <stdio.h>
#include <queue>

//RNAMake Headers
#include "motif/motif_state.h"
#include "motif_search/motif_state_search_node.fwd.h"

namespace motif_search {

class MotifStateSearchNode {
public:
    inline
    MotifStateSearchNode(
            motif::MotifStateOP const & ref_state,
            MotifStateSearchNodeOP const & parent,
            int parent_end_index,
            int ntype) :
            ref_state_(ref_state),
            parent_(parent),
            parent_end_index_(parent_end_index),
            ntype_(ntype),
            node_type_usages_(Ints()) {
        ss_score_ = ref_state_->score();
        size_ = ref_state_->size();
        if (parent_ == nullptr) { level_ = 1; }
        else {
            level_ = parent->level_;
            ss_score_ += parent_->ss_score_;
            size_ += parent_->size_;
        }
        cur_state_ = std::make_shared<motif::MotifState>(*ref_state_);
        score_ = 1000;
    }

    inline
    MotifStateSearchNode(
            MotifStateSearchNode const & n) :
            ref_state_(n.ref_state_),
            parent_(n.parent_),
            parent_end_index_(n.parent_end_index_),
            ntype_(n.ntype_),
            node_type_usages_(n.node_type_usages_),
            cur_state_(std::make_shared<motif::MotifState>(*n.cur_state_)),
            ss_score_(n.ss_score_),
            level_(n.level_),
            size_(n.size_),
            score_(n.score_),
            center_(n.center_) {}

public:
    MotifStateSearchNode
    copy() {
        MotifStateSearchNode new_n(ref_state_, parent_, parent_end_index_, ntype_);
        new_n.cur_state_ = std::make_shared<motif::MotifState>(*cur_state_);
        new_n.score_ = score_;
        new_n.ss_score_ = ss_score_;
        new_n.level_ = level_;
        new_n.node_type_usages_ = node_type_usages_;
        return new_n;
    }

    inline
    void
    calc_center() {
        center_.x(0);
        center_.y(0);
        center_.z(0);
        for (auto const & b : cur_state_->beads()) {
            center_ += b;
        }
        center_ /= cur_state_->beads().size();
    }

    inline
    int
    node_type_usage(
            int i) {
        if (i == -1) { return 0; }
        return node_type_usages_[i];
    }

    inline
    void
    update(int parent_end_index) {
        if (parent_ == nullptr) { return; }
        level_ = parent_->level_ + 1;
        parent_end_index_ = parent_end_index;
        ss_score_ = ref_state_->score() + parent_->ss_score_;
        size_ = ref_state_->size() + parent_->size_;
        if (ntype_ == -1) { return; }
        node_type_usages_ = parent_->node_type_usages_;
        node_type_usages_[ntype_] += 1;

    }


    inline
    void
    replace_ms(
            motif::MotifStateOP const & ms,
            int ntype) {
        ref_state_ = ms;
        cur_state_ = ms;
        ntype_ = ntype;

    }

    inline
    void
    setup_node_type_usage(
            int size) {

        node_type_usages_ = Ints(size);

    }

public: //getters

    inline
    MotifStateSearchNodeOP const &
    parent() { return parent_; }

    inline
    motif::MotifStateOP &
    cur_state() { return cur_state_; }

    inline
    motif::MotifStateOP const &
    ref_state() { return ref_state_; }

    inline
    int
    level() { return level_; }

    inline
    int
    size() { return size_; }

    inline
    int
    ntype() { return ntype_; }

    inline
    float
    ss_score() { return ss_score_; }

    inline
    float
    score() { return score_; }

    inline
    int
    parent_end_index() { return parent_end_index_; }

    inline
    math::Point const &
    center() { return center_; }


public: //setters

    inline
    void
    score(
            float score) { score_ = score; }

    inline
    void
    level(
            int level) { level_ = level; }

    inline
    void
    parent(
            MotifStateSearchNodeOP const & parent) { parent_ = parent; }

private:
    motif::MotifStateOP ref_state_, cur_state_;
    MotifStateSearchNodeOP parent_;
    math::Point center_;
    Ints node_type_usages_;
    int parent_end_index_, level_, size_, ntype_;
    float ss_score_, score_;

};


struct MotifStateSearchNodeOPCompare {
    bool
    operator()(
            MotifStateSearchNodeOP const & node1,
            MotifStateSearchNodeOP const & node2) {

        if (node1->score() > node2->score()) { return true; }
        else { return false; }
    }
};


typedef std::priority_queue<MotifStateSearchNodeOP, MotifStateSearchNodeOPs, MotifStateSearchNodeOPCompare> MotifStateSearchNodeQueue;

}

#endif /* defined(__RNAMake__motif_state_search_node__) */
