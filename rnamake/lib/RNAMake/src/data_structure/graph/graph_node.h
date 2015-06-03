//
//  graph_node.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__graph_node__
#define __RNAMake__graph_node__

#include <stdio.h>

#include "base/types.h"
#include "data_structure/graph/graph_node.fwd.h"

enum class GraphNodeType { GraphNodeTypeStatic, GraphNodeTypeDynamic };

template <typename DataType>
class GraphNode {
public:
    GraphNode(
        DataType const & data,
        int index,
        int level,
        GraphNodeType type = GraphNodeType::GraphNodeTypeDynamic,
        size_t n_connections = 0):
    data_(data),
    level_(level),
    index_(index),
    type_(type),
    connections_(GraphConnectionOPs<DataType>(n_connections))
    {}
    
public:
    
    inline
    void
    add_connection(
        GraphConnectionOP<DataType> const & connection,
        int pos=-1) {
        
        if(pos == -1 && type_ == GraphNodeType::GraphNodeTypeStatic) {
            throw std::runtime_error("attempted to resize children array in NodeTypeStatic GraphNode");
        }
        
        if(pos == -1) {
            connections_.push_back(connection);
            return;
        }
        
        if(pos >= connections_.size()) {
            throw std::runtime_error("cannot add child at position");
        }
        
        if(connections_[pos] != nullptr) {
            throw std::runtime_error("attempted to add child in a position that is already full");
        }
        
        connections_[pos] = connection;
    }
    
    inline
    Ints
    available_children_pos() {
        Ints pos;
        int i = 0;
        for(auto const & c : connections_) {
            if(c == nullptr) { pos.push_back(i); }
        }
        return pos;
    }

    
public:
    inline
    int
    index() const { return index_; }
    
    inline
    int
    level() const { return level_; }

    
private:
    DataType data_;
    GraphConnectionOPs<DataType> connections_;
    GraphNodeType type_;
    int level_, index_;
    
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
    
private:
    GraphNodeOP<DataType> node_1_, node_2_;
    int end_index_1_, end_index_2_;
};


#endif /* defined(__RNAMake__graph_node__) */
