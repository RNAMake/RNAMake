//
//  tree_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__tree_node__
#define __RNAMake__tree_node__

#include <stdio.h>
#include <vector>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <algorithm>

//RNAMake Headers
#include "base/types.hpp"
#include "data_structure/tree/tree_node.fwd.hh"

namespace data_structure {
namespace tree {

class TreeException : public std::runtime_error {
public:
    TreeException(
            String const & message) :
            std::runtime_error(message) {}

};

template<typename DataType>
class TreeNode {

public:
    TreeNode(
            DataType const & data,
            int index,
            int level,
            size_t n_children = 0) :
            data_(data),
            level_(level),
            index_(index),
            children_(TreeNodeOPs<DataType>(n_children)) {}

public: //Vitrual functions need to be implemented in derived clases

    virtual
    void
    add_child(
            TreeNodeOP <DataType> const &,
            int pos = -1) = 0;

    virtual
    void
    remove_child(
            TreeNodeOP <DataType> const &) = 0;


    virtual
    bool
    leaf() = 0;


public:
    inline
    Indexes
    available_children_pos() const {
        Indexes pos;
        int i = 0;
        for (auto const & c : children_) {
            if (c == nullptr) { pos.push_back(i); }
            i++;
        }
        return pos;
    }

    inline
    int
    available_pos(int pos) {
        if (children_.size() <= pos) { return 0; }
        if (children_[pos] != nullptr) { return 0; }
        return 1;
    }

    inline
    void
    unset_children() {
        for (auto & c : children_) { c = nullptr; }
    }


public: //getters
    inline
    int
    index() const { return index_; }

    inline
    int
    level() const { return level_; }

    inline
    DataType const &
    data() const { return data_; }

    inline
    DataType &
    data() { return data_; }

    inline
    TreeNodeOPs <DataType> const &
    children() { return children_; }

    inline
    TreeNodeOP <DataType> const &
    parent() { return parent_; }

    inline
    int
    parent_index() {
        if (parent_ == nullptr) {
            return -1;
        }
        return parent_->index();
    }

    inline
    int
    parent_end_index() {
        if (parent_ == nullptr) { return -1; }
        int i = -1;
        for (auto const & c : parent_->children()) {
            i++;
            if (c == nullptr) { continue; }
            if (c->index() == index_) { return i; }
        }
        return -1;

    }

public: //setters

    inline
    void
    parent(
            TreeNodeOP <DataType> const & p) {
        parent_ = p;
    }

    inline
    void
    index(
            int index) {
        index_ = index;
    }

protected:
    DataType data_;
    TreeNodeOPs <DataType> children_;
    TreeNodeOP <DataType> parent_;
    int level_, index_;

};

template<typename DataType>
class TreeNodeDynamic : public TreeNode<DataType> {
public:
    inline
    TreeNodeDynamic(
            DataType const & data,
            int index,
            int level) :
            TreeNode<DataType>(data, index, level, 0) {}

public:

    inline
    void
    add_child(
            TreeNodeOP <DataType> const & c,
            int pos = -1) {
        this->children_.push_back(c);
    }

    inline
    void
    remove_child(
            TreeNodeOP <DataType> const & child) {

        int found = 0;
        for (auto & c : this->children_) {
            if (c == child) {
                found = 1;
                break;
            }
        }

        if (!found) {
            throw TreeException("tried to remove child but is not present in node");
        }

        this->children_.erase(std::remove(this->children_.begin(), this->children_.end(),
                                          child), this->children_.end());

    }

    inline
    bool
    leaf() {
        if (this->children_.size() == 0) { return true; }
        else { return false; }
    }

};

template<typename DataType>
class TreeNodeStatic : public TreeNode<DataType> {
public:
    inline
    TreeNodeStatic(
            DataType const & data,
            int index,
            int level,
            int n_children) :
            TreeNode<DataType>(data, index, level, n_children) {}

public:

    inline
    void
    add_child(
            TreeNodeOP <DataType> const & child,
            int pos = -1) {

        if (pos == -1) {
            throw TreeException("attempted to resize children array in TreeNodeStatic");
        }

        if (pos >= this->children_.size()) {
            throw TreeException("cannot add child at position");
        }

        if (this->children_[pos] != nullptr) {
            throw TreeException("attempted to add child in a position that is already full");
        }

        this->children_[pos] = child;
    }

    inline
    void
    remove_child(
            TreeNodeOP <DataType> const & child) {

        for (auto & c : this->children_) {
            if (c == child) {
                c = nullptr;
                return;
            }
        }

        throw TreeException("tried to remove child but is not present in node");

    }

    inline
    bool
    leaf() {
        for (auto const & c : this->children_) {
            if (c != nullptr) { return false; }
        }

        return true;
    }
};

}
}









#endif /* defined(__RNAMake__tree_node__) */
