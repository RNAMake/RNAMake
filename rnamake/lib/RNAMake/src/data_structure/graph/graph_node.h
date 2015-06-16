
//  graph_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__graph_node__
#define __RNAMake__graph_node__

#include <stdio.h>
#include <iostream>

#include "base/types.h"
#include "data_structure/graph/graph_node.fwd.h"

class GraphException : public std::runtime_error {
public:
    GraphException(
        String const & message) :
    std::runtime_error("Graph Exception: " + message)
    {}
    
};


template <typename DataType>
class GraphNode {

public:
    GraphNode(
        DataType const & data,
        int index,
        int level,
        size_t n_connections = 0):
    data_(data),
    level_(level),
    index_(index),
    connections_(GraphConnectionOPs<DataType>(n_connections))
    {}
    
public:
    
    virtual
    inline
    void
    add_connection(
        GraphConnectionOP<DataType> const &,
        int pos=-1) = 0;
    
    virtual
    inline
    void
    remove_connection(
        GraphConnectionOP<DataType> const & connection) = 0;
    
    inline
    Ints
    available_children_pos() const {
        Ints pos;
        int i = 0;
        for(auto const & c : connections_) {
            if(c == nullptr) { pos.push_back(i); }
            i++;
        }
        return pos;
    }
    
    inline
    int
    available_pos(int pos) {
        if(connections_.size() <= pos) { return 0; }
        if(connections_[pos] != nullptr) { return 0; }
        return 1;
    }

    
public:
    inline
    GraphNodeOP<DataType> const &
    parent() {
        for(auto const & c : connections_) {
            if(c->partner(index_)->index_ < index_) { return c->partner(index_); }
        }
        
        throw GraphException("could not find parent for node");
    }
    
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
    GraphConnectionOPs<DataType> const &
    connections() { return connections_; }
    
    
protected:
    DataType data_;
    GraphConnectionOPs<DataType> connections_;
    int level_, index_;
    
};


template <typename DataType>
class GraphNodeDynamic : public GraphNode<DataType> {
public:
    inline
    GraphNodeDynamic(
        DataType const & data,
        int index,
        int level):
    GraphNode<DataType>(data, index, level, 0) {}
    
public:
    
    inline
    void
    add_connection(
        GraphConnectionOP<DataType> const & connection,
        int pos=-1) {
        
        this->connections_.push_back(connection);
    }
    
    inline
    void
    remove_connection(
        GraphConnectionOP<DataType> const & connection)  {
        
        int found = 0;
        for(auto & c : this->connections_) {
            if(c == connection) {
                found = 1;
                break;
            }
        }
        
        if(!found) {
            throw GraphException("tried to remove connection but is not present in node");
        }
        
        this->connections_.erase(std::remove(this->connections_.begin(), this->connections_.end(),
                                connection), this->connections_.end());

    }

    
};


template <typename DataType>
class GraphNodeStatic  : public GraphNode<DataType> {
public:
    inline
    GraphNodeStatic(
        DataType const & data,
        int index,
        int level,
        int n_children):
    GraphNode<DataType>(data, index, level, n_children) {}
    
public:
    
    inline
    void
    add_connection(
        GraphConnectionOP<DataType> const & connection,
        int pos=-1) {
                
        if(pos == -1) {
            throw GraphException("attempted to resize children array in NodeTypeStatic GraphNode");
        }
        
        if(pos >= this->connections_.size()) {
            throw GraphException("cannot add child at position");
        }
        
        if(this->connections_[pos] != nullptr) {
            throw GraphException("attempted to add child in a position that is already full");
        }
        
        this->connections_[pos] = connection;
    }
    
    inline
    void
    remove_connection(
        GraphConnectionOP<DataType> const & connection)  {
        
        for(auto & c : this->connections_) {
            if(c == connection) {
                c = nullptr;
                return;
            }
        }
        
        throw GraphException("tried to remove connection but is not present in node");
        
    }

    
    
};


template <typename DataType>
class GraphConnection {
public:
    GraphConnection(
        GraphNodeOP<DataType> node_1,
        GraphNodeOP<DataType> node_2,
        int end_index_1,
        int end_index_2):
    node_1_(node_1),
    node_2_(node_2),
    end_index_1_(end_index_1),
    end_index_2_(end_index_2)
    {}
    
public:
    void
    disconnect() {
        node_1_ = nullptr;
        node_2_ = nullptr;
    }
    
    GraphNodeOP<DataType> const &
    partner(
        int i) {
        if     (i == node_1_->index()) { return node_2_; }
        else if(i == node_2_->index()) { return node_1_; }
        else { throw GraphException("cannot call partner with node not in connection"); }
    }
    
    
private:
    GraphNodeOP<DataType> node_1_, node_2_;
    int end_index_1_, end_index_2_;
};


#endif /* defined(__RNAMake__graph_node__) */