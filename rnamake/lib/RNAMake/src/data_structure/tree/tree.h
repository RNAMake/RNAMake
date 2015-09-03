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

//RNAMAke Headers
#include "base/types.h"
#include "data_structure/tree/tree_node.fwd.hh"
#include "data_structure/tree/tree_node.h"


template <typename DataType>
class Tree {
public:
    
    Tree():
    index_(0),
    level_(0)
    {}
    
    virtual
    ~Tree() {
        for(int i = 0; i < nodes_.size(); i++){
            nodes_[i]->unset_children();
        }
    }
    
public:
    typedef typename TreeNodeOPs<DataType>::iterator iterator;
    typedef typename TreeNodeOPs<DataType>::const_iterator const_iterator;
    
    iterator begin() { return nodes_.begin(); }
    iterator end()   { return nodes_.end(); }
    
    const_iterator begin() const { return nodes_.begin(); }
    const_iterator end()   const { return nodes_.end(); }
    
public:
    
    
    inline
    TreeNodeOP<DataType> const &
    get_node(
        int index) {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw TreeException("cannot find node with index");
    }
    
    inline
    void
    remove_node(TreeNodeOP<DataType> const & n) {
        if(n->parent() != nullptr) { n->parent()->remove_child(n); }
        n->parent(nullptr);
        n->unset_children();
        nodes_.erase(std::remove(nodes_.begin(), nodes_.end(), n), nodes_.end());
        last_node_ = nodes_.back();
    }

    inline
    void
    remove_node(int pos) { return remove_node(get_node(pos)); }
    
public: //getters
    
    inline
    size_t
    size() { return nodes_.size(); }
    
    inline
    TreeNodeOP<DataType> const &
    last_node() { return last_node_; }
    
    
protected:
    TreeNodeOPs<DataType> nodes_;
    TreeNodeOP<DataType> last_node_;
    int level_, index_;
    
};

template <typename DataType>
class TreeDynamic : public Tree<DataType> {
public:
    TreeDynamic(): Tree<DataType>() {}
    
    ~TreeDynamic() {
        for(int i = 0; i < this->nodes_.size(); i++){
            this->nodes_[i]->unset_children();
        }
    }
    
public:
   
    inline
    int
    add_data(
        DataType const & data,
        int parent_index=-1) {
        
        auto parent = this->last_node_;
        if(parent_index != -1) { parent = this->get_node(parent_index); }
        auto n = std::make_shared<TreeNodeDynamic<DataType>>(data, this->index_, this->level_);
        
        if(parent != nullptr) {
            parent->add_child(n);
            n->parent(parent);
        }
        
        this->nodes_.push_back(n);
        this->index_++;
        this->last_node_ = n;
        return this->index_-1;
    
    }
    
};


template <typename DataType>
class TreeStatic : public Tree<DataType> {
public:
    TreeStatic(): Tree<DataType>() {}
    
    ~TreeStatic() {
        for(int i = 0; i < this->nodes_.size(); i++){
            this->nodes_[i]->unset_children();
        }
    }
    
public:
    
    inline
    int
    add_data(
        DataType const & data,
        int n_children = 1,
        int parent_index=-1,
        int parent_child_index=-1) {
        
        auto parent = this->last_node_;
        if(parent_index != -1) { parent = this->get_node(parent_index); }
        auto n = std::make_shared<TreeNodeStatic<DataType>>(data, this->index_, this->level_, n_children);
        
        if(parent != nullptr) {
            if(parent_child_index == -1) {
                auto all_pos = parent->available_children_pos();
                if(all_pos.size() == 0) {
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
        return this->index_-1;
    }
    
    inline
    Ints
    get_available_pos(
        TreeNodeOP<DataType> const & n,
        int pos=-1) {
        
        if(pos == -1) {
            return n->available_children_pos();
        }
        else {
            if(n->available_pos(pos) == 0) {
                throw TreeException("tree pos is not available");
            }
            return Ints { pos } ;
        }
        
    }
    
    
};


#endif /* defined(__RNAMake__tree__) */

























