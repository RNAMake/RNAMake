//
//  graph.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__graph__
#define __RNAMake__graph__

#include <stdio.h>

#include "data_structure/graph/graph_node.h"
#include "data_structure/graph/graph_node.fwd.h"

template <typename DataType>
class GraphIterator;


template <typename DataType>
class Graph {
public:
    Graph():
    index_(0),
    level_(0) {}
    
    ~Graph() {}
    
public:
    
    typedef GraphIterator<DataType> iterator;
    friend class GraphIterator<DataType>;
    iterator begin() const;
    iterator end() const;
    
public:
    
    inline
    size_t
    size() { return nodes_.size(); }
    
    inline
    GraphNodeOP<DataType> const &
    get_node(
        int index) {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw GraphException("cannot find node with index");
    }

    inline
    GraphNodeOPs<DataType> const &
    nodes() {
        return nodes_;
    }
    
    
protected:
    GraphNodeOP<DataType> last_node_;
    GraphNodeOPs<DataType> nodes_;
    GraphConnectionOPs<DataType> connections_;
    int index_, level_;
};

template <typename DataType>
class GraphDynamic : public Graph<DataType> {
public:
    GraphDynamic(): Graph<DataType>() {}
   
public:
    
    inline
    int
    add_data(
        DataType const & data,
        int parent_index = -1) {
        
        GraphNodeOP<DataType> parent = this->last_node_;
        if(parent_index != -1) { parent = this->get_node(parent_index);}
        
        auto n = std::make_shared<GraphNodeDynamic<DataType>>(data, this->index_, this->level_);
        
        if(parent != nullptr) {
            auto c = std::make_shared<GraphConnection<DataType>>(parent, n, 0, 0);
            parent->add_connection(c);
            n->add_connection(c);
        }
        
    
        this->nodes_.push_back(n);
        this->index_++;
        this->last_node_ = n;
        return this->index_-1;
    }
    
    
};

template <typename DataType>
class GraphStatic : public Graph<DataType> {
public:
    GraphStatic(): Graph<DataType>() {}
    
public:
    
    inline
    int
    add_data(
        DataType const & data,
        int parent_index = -1,
        int parent_pos = -1,
        int child_pos = -1,
        int n_children = 0) {
        
        GraphNodeOP<DataType> parent = this->last_node_;
        auto n = std::make_shared<GraphNodeStatic<DataType>>(data, this->index_, this->level_,
                                                             n_children);

        if(parent_index != -1) { parent = this->get_node(parent_index); }
        if(parent != nullptr) {
            check_pos_is_valid(parent, parent_pos);
            check_pos_is_valid(n, child_pos);
            auto c = std::make_shared<GraphConnection<DataType>>(parent, n, parent_pos, child_pos);
            parent->add_connection(c, parent_pos);
            n->add_connection(c, child_pos);
        }
        
        this->nodes_.push_back(n);
        this->index_++;
        this->last_node_ = n;
        return this->index_-1;
        
        
    }
    
    inline
    void
    check_pos_is_valid(
        GraphNodeOP<DataType> const & n,
        int & pos) {
        
        if(pos == -1) {
            Ints avail_pos = n->available_children_pos();
            if(avail_pos.size() == 0) {
                throw GraphException("cannot add connection to node, has not available ends");
            }
            pos = avail_pos[0];
            return;
        }
        
        else {
            if(n->available_pos(pos) == 0) {
                throw GraphException("graph pos is not availabe");
            }
        }
        
    }
    
};

template <typename DataType>
class GraphIterator {
public:
    friend class Graph<DataType>;
    
    GraphIterator() {}
    
public:
    GraphIterator & operator++ ();
    GraphIterator operator++ (int);
    
    bool operator== (const GraphIterator& rhs) const;
    bool operator!= (const GraphIterator& rhs) const;
    GraphNodeOP<DataType>& operator* ();
    
private:
    GraphNodeOP<DataType> node_ptr_;
    Graph<DataType>* graph_;
    
    GraphIterator(
        GraphNodeOP<DataType> const & node,
        Graph<DataType>* graph):
    node_ptr_(node),
    graph_(graph)
    {}
    
};

template <typename DataType>
typename Graph<DataType>::iterator
Graph<DataType>::begin() const {
    return iterator(nodes_[0], this);
}

template <typename DataType>
typename Graph<DataType>::iterator
Graph<DataType>::end() const {
    return iterator(nullptr, nullptr);
}

template <typename DataType>
GraphNodeOP<DataType> &
GraphIterator<DataType>::operator* () {
    return node_ptr_;
}

template <typename DataType>
GraphIterator<DataType> &
GraphIterator<DataType>::operator++() {
    if(node_ptr_->index()+1 == (int)graph_->size()) {
        node_ptr_ = nullptr;
    }
    else{
        node_ptr_ = graph_->get_node(node_ptr_->index()+1);
    }
    
    return *this;
}

template <typename DataType>
bool
GraphIterator<DataType>::operator== (
    GraphIterator<DataType> const& rhs) const {
    return node_ptr_ == rhs.node_ptr_;
}

template <typename DataType>
bool
GraphIterator<DataType>::operator!= (
    GraphIterator<DataType> const& rhs) const {
    return node_ptr_ != rhs.node_ptr_;
}




#endif /* defined(__RNAMake__graph__) */
