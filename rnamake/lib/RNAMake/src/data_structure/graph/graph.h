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
class Graph {
public:
    Graph():
    index_(0),
    level_(0) {}
    
    ~Graph() {}
    
public:
    
    
    inline
    GraphNodeOP<DataType> const &
    get_node(
        int index) {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw GraphException("cannot find node with index");
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





#endif /* defined(__RNAMake__graph__) */
