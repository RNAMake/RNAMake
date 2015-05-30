//
//  node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__node__
#define __RNAMake__node__

#include <stdio.h>
#include <vector>
#include <memory>
#include <iostream>

//RNAMake Headers
#include "base/types.h"
#include "data_structure/tree/node.fwd.h"

enum NodeChildType { NodeTypeStatic, NodeTypeDynamic };


template <class DataType>
class Node {
public:
    
    inline
    Node(
        DataType const & data,
        int index,
        int level,
        NodeChildType type = NodeTypeDynamic,
        size_t const n_children=0):
    index_(index),
    level_(level),
    children_(std::vector<Node<DataType>*>(n_children)),
    parent_(nullptr),
    type_(type),
    data_(data)
    {}
    
    virtual ~Node() {}
    
public:
    
    inline
    void
    add_child(
        Node<DataType>* const & child,
        int pos=-1) {
        
        if(pos == -1 && type_ == NodeTypeStatic) {
            throw std::runtime_error("attempted to resize children array in NodeTypeStatic Node");
        }
        
        if(pos == -1) {
            children_.push_back(child);
            return;
        }
        
        if(pos >= children_.size()) {
            throw std::runtime_error("cannot add child at position");
        }
        
        if(children_[pos] != nullptr) {
            throw std::runtime_error("attempted to add child in a position that is already full");
        }
        
        children_[pos] = child;
    }
   
    inline
    void
    unset_child(
        Node<DataType>* const & child) {
        
        for(auto & c : children_) {
            if(c == child) {
                c = nullptr;
                return;
            }
        }
        
        throw std::runtime_error("cannot remove this child it is not a child of this node");
    }
    
    
    inline
    Ints
    available_children_pos() {
        Ints pos;
        int i = 0;
        for(auto const & c : children_) {
            if(c == nullptr) { pos.push_back(i); }
        }
        return pos;
    }

    inline
    int
    parent_index() {
        if(parent_ == nullptr) {
            return -1;
        }
        return parent_->index_;
    }
    

public:
    
    virtual
    inline
    void
    set_parent(
        Node<DataType>* const & nparent) {
        parent_ = nparent;
    }
    
public: //getters
    
    inline
    Node<DataType>* const &
    parent() { return parent_; }
    
    inline
    DataType const &
    data() { return data_; }
    
    inline
    std::vector<Node<DataType>*> const &
    children() {
        return children_;
    }
    
    inline
    int
    index() { return index_; }
    
    inline
    int
    level() { return level_; }

    
    
    
protected:
    int index_, level_;
    Node<DataType>* parent_;
    std::vector<Node<DataType>*> children_;
    DataType data_;
    NodeChildType type_;

};





#endif /* defined(__RNAMake__node__) */
