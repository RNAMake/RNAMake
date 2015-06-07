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
#include "data_structure/tree/node.fwd.h"
#include "data_structure/tree/node.h"


template <typename DataType>
class TreeIterator;

template <typename DataType>
class Tree {
public:
    
    Tree(
         NodeChildType const & node_type):
    nodes_ (NodeOPs<DataType>() ),
    last_node_(nullptr),
    level_(0),
    index_(0),
    node_type_(node_type){
    }
    
    Tree():
    nodes_ (NodeOPs<DataType>()),
    last_node_(nullptr),
    level_(0),
    index_(0),
    node_type_(NodeTypeDynamic) {
    }
    
    
    ~Tree() {
        for(int i = 0; i < nodes_.size(); i++){
            nodes_[i]->unset_children();
        }
    }
    
public:
    typedef TreeIterator<DataType> iterator;
    friend class TreeIterator<DataType>;
    iterator begin() const;
    iterator end() const;
    
public:
    
    inline
    int
    add_data(
        DataType const & data,
        int n_children=0,
        int parent_index=-1,
        int parent_child_pos=-1) {
        auto n = std::make_shared<Node<DataType>>(data, index_, level_, node_type_, n_children);
        NodeOP<DataType> parent = last_node_;
        if(parent_index != -1) { parent = get_node(parent_index);}
        if(parent == nullptr) {
            nodes_.push_back(n);
            index_++;
            last_node_ = n;
            return 0;
        }
        else if (parent_child_pos == -1 && node_type_ == NodeTypeStatic) {
            Ints avail_child_pos = parent->available_children_pos();
            if(avail_child_pos.size() == 0) {
                throw std::runtime_error("cannot use node as parent as it has not spots for children\n");
            }
            parent_child_pos = avail_child_pos[0];
        }
     
        parent->add_child(n, parent_child_pos);
        n->set_parent(parent);
        
        nodes_.push_back(n);
        index_++;
        last_node_ = n;
    
        return n->index();
        
    }
    
    inline
    NodeOP<DataType> const &
    get_node(
        int index) {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw std::runtime_error("cannot find node with index");
    }
    
    inline
    int
    get_index(
        DataType const & data) {
        int index = -1;
        int found = 0;
        for(auto const & n : nodes_) {
            if(n->data() == data) {
                index = n->index();
                found++;
            }
        }
        
        if(index == -1) {
            throw std::runtime_error("cannot find index of given data in Tree");
        }
        if(found > 1) {
            throw std::runtime_error("cannot get index since there is duplicate data in Tree");
        }
        
        return index;
    }
    
    inline
    DataType const &
    get_data(
        int index) {
        return get_node(index)->data();
    }
    
    inline
    void
    remove_node(int pos) {
        auto node = get_node(pos);
        for(auto & c : node->children()) {
            if(c == nullptr) { continue; }
            c->set_parent(nullptr);
            remove_node(c->index());
            node->unset_child(c);
        }
        
       
        last_node_ = node->parent();
        node->set_parent(nullptr);
        nodes_.erase(std::remove(nodes_.begin(), nodes_.end(), node), nodes_.end());
    
    }
    
public: //getters
    
    inline
    size_t
    size() { return nodes_.size(); }
    
    inline
    NodeOPs<DataType> const &
    nodes() const {
        return nodes_;
    }
    
    
protected:
    NodeOPs<DataType> nodes_;
    NodeOP<DataType> last_node_;
    int level_, index_;
    NodeChildType node_type_;
    
};

template <typename DataType>
class TreeIterator {
public:
    friend class Tree<DataType>;
    
    TreeIterator() {}
    
public:
    TreeIterator & operator++ ();
    TreeIterator operator++ (int);
    
    bool operator== (const TreeIterator& rhs) const;
    bool operator!= (const TreeIterator& rhs) const;
    NodeOP<DataType>& operator* ();
    
private:
    NodeOP<DataType>node_ptr_;
    std::queue<NodeOP<DataType>> queue_;
    
    TreeIterator(
        NodeOP<DataType> const & node):
        node_ptr_(node),
        queue_(std::queue<NodeOP<DataType>>())
    {}
};


template <typename DataType>
typename Tree<DataType>::iterator
Tree<DataType>::begin() const {
    return iterator(nodes_[0]);
}

template <typename DataType>
typename Tree<DataType>::iterator
Tree<DataType>::end() const {
    return iterator(nullptr);
}

template <typename DataType>
NodeOP<DataType>&
TreeIterator<DataType>::operator* () {
    return node_ptr_;
}

template <typename DataType>
TreeIterator<DataType>&
TreeIterator<DataType>::operator++() {
    auto children = node_ptr_->children();
    for(auto const & c : children) {
        if(c != nullptr) { queue_.push(c); }
    }
    
    if(!queue_.empty()) {
        node_ptr_ = queue_.front();
        queue_.pop();
    }
    else {
        node_ptr_ = nullptr;
    }
    
    
    return *this;
}

template <typename DataType>
bool
TreeIterator<DataType>::operator== (
    TreeIterator<DataType> const& rhs) const {
    return node_ptr_ == rhs.node_ptr_;
}

template <typename DataType>
bool
TreeIterator<DataType>::operator!= (
    TreeIterator<DataType> const& rhs) const {
    return node_ptr_ != rhs.node_ptr_;
}




#endif /* defined(__RNAMake__tree__) */





























