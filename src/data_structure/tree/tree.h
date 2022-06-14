//
//  graph.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__tree__
#define __RNAMake__tree__

#include <stdio.h>
#include <memory>
#include <queue>
#include <algorithm>
#include <cassert>

//RNAMAke Headers
#include "base/types.hpp"
#include "data_structure/tree/tree_node.fwd.hh"
#include "data_structure/tree/tree_node.h"

namespace data_structure {
namespace tree {

template<typename DataType>
class Tree {
public:

    Tree() :
            index_(0),
            level_(0) {}

    virtual
    ~Tree() {
        for (int i = 0; i < nodes_.size(); i++) {
            nodes_[i]->unset_children();
        }
    }

public:
    typedef typename TreeNodeOPs<DataType>::iterator iterator;
    typedef typename TreeNodeOPs<DataType>::const_iterator const_iterator;

    iterator begin() { return nodes_.begin(); }

    iterator end() { return nodes_.end(); }

    const_iterator begin() const { return nodes_.begin(); }

    const_iterator end() const { return nodes_.end(); }

public:


    inline
    TreeNodeOP<DataType> const &
    get_node(
            int index) {

        for (auto const & n : nodes_) {
            if (n->index() == index) { return n; }
        }

        throw TreeException("cannot find node with index: " + std::to_string(index));
    }

    inline
    void
    remove_node(TreeNodeOP<DataType> const & n) {
        if (n->parent() != nullptr) { n->parent()->remove_child(n); }
        n->parent(nullptr);
        n->unset_children();
        nodes_.erase(std::remove(nodes_.begin(), nodes_.end(), n), nodes_.end());
        if (nodes_.size() > 0) {
            last_node_ = nodes_.back();
        } else {
            last_node_ = nullptr;
        }
    }

    inline
    void
    remove_node(int pos) { return remove_node(get_node(pos)); }

    inline
    void
    remove_level(
            int level) {

        int pos = 0;
        int removed = 1;
        while (removed) {
            removed = 0;
            for (auto const & n : nodes_) {
                if (n->level() >= level && n->leaf()) {
                    remove_node(n);
                    removed = 1;
                    break;
                }
            }
        }

    }

    inline
    void
    increase_level() { level_ += 1; }

    inline
    void
    decrease_level() {
        level_ -= 1;
        assert(level_ > -1 && "level has to be positive");
    }

public: //getters

    inline
    size_t
    size() { return nodes_.size(); }

    inline
    TreeNodeOP<DataType> const &
    last_node() { return last_node_; }

    inline
    int
    level() { return level_; }


protected:
    TreeNodeOPs<DataType> nodes_;
    TreeNodeOP<DataType> last_node_;
    int level_, index_;

};

template<typename DataType>
class TreeDynamic : public Tree<DataType> {
public:
    TreeDynamic() : Tree<DataType>() {}

    ~TreeDynamic() {
        for (int i = 0; i < this->nodes_.size(); i++) {
            this->nodes_[i]->unset_children();
        }
    }

public:

    inline
    int
    add_data(
            DataType const & data,
            int parent_index = -1) {

        auto parent = this->last_node_;
        if (parent_index != -1) { parent = this->get_node(parent_index); }
        auto n = std::make_shared<TreeNodeDynamic<DataType>>(data, this->index_, this->level_);

        if (parent != nullptr) {
            parent->add_child(n);
            n->parent(parent);
        }

        this->nodes_.push_back(n);
        this->index_++;
        this->last_node_ = n;
        return this->index_ - 1;

    }

};


template<typename DataType>
class TreeStatic : public Tree<DataType> {
public:
    TreeStatic() : Tree<DataType>() {}


    TreeStatic(
            TreeStatic<DataType> const & t) {

        for (auto const & n : t) {
            auto n_new = std::make_shared<TreeNodeStatic<DataType>>(n->data(), n->index(),
                                                                    n->level(), n->children().size());
            if (n->parent() != nullptr) {
                auto parent = this->get_node(n->parent_index());
                parent->add_child(n_new, n->parent_end_index());
                n_new->parent(parent);
            }

            this->nodes_.push_back(n_new);

        }

        if (t.last_node_ != nullptr) {
            this->last_node_ = this->get_node(t.last_node_->index());
        }

        this->level_ = t.level_;
        this->index_ = t.index_;

    }

    ~TreeStatic() {
        for (int i = 0; i < this->nodes_.size(); i++) {
            this->nodes_[i]->unset_children();
        }
    }

public:

    inline
    int
    add_data(
            DataType const & data,
            int n_children = 1,
            int parent_index = -1,
            int parent_child_index = -1) {

        auto parent = this->last_node_;
        if (parent_index != -1) { parent = this->get_node(parent_index); }
        auto n = std::make_shared<TreeNodeStatic<DataType>>(data, this->index_, this->level_, n_children);

        if (parent != nullptr) {
            if (parent_child_index == -1) {
                auto all_pos = parent->available_children_pos();
                if (all_pos.size() == 0) {
                    throw TreeException("cannot add child to this parent it has no open spots");
                }
                parent_child_index = all_pos[0];
            }
            parent->add_child(n, parent_child_index);
            n->parent(parent);
        }

        this->nodes_.push_back(n);
        this->index_++;
        this->last_node_ = n;
        return this->index_ - 1;
    }

    inline
    Indexes
    get_available_pos(
            TreeNodeOP<DataType> const & n,
            int pos = -1) {

        if (pos == -1) {
            return n->available_children_pos();
        } else {
            if (n->available_pos(pos) == 0) {
                throw TreeException("tree pos is not available");
            }
            return Indexes {pos};
        }

    }

public:
    inline
    void
    index(
            int n_index) {
        this->index_ = n_index;
    }
};

}
}

#endif /* defined(__RNAMake__tree__) */

























