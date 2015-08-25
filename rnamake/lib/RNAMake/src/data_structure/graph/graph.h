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
#include <sstream> 
#include <typeinfo>

#include "data_structure/graph/graph_node.h"
#include "data_structure/graph/graph_node.fwd.h"


template <typename DataType>
class Graph {
public:
    Graph():
    index_(0),
    level_(0) {}
    
    virtual
    ~Graph() {
        for(int i = 0; i < nodes_.size(); i++){
            nodes_[i]->unset_connections();
        }
    }
    
public:
    
    typedef typename GraphNodeOPs<DataType>::iterator iterator;
    typedef typename GraphNodeOPs<DataType>::const_iterator const_iterator;

    iterator begin() { return nodes_.begin(); }
    iterator end()   { return nodes_.end(); }
    
    const_iterator begin() const { return nodes_.begin(); }
    const_iterator end()   const { return nodes_.end(); }
        
public:
    
    inline
    size_t
    size() const { return nodes_.size(); }
    
    inline
    GraphNodeOP<DataType> const &
    get_node(
        int index) const {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw GraphException("cannot find node with index");
    }


    inline
    void
    increase_level() { level_ += 1; }
    
    inline
    void
    decrease_level() { level_ -= 1; }
    
    
public: //getters
    
    inline
    GraphNodeOPs<DataType> const &
    nodes() {
        return nodes_;
    }
    
    inline
    GraphNodeOP<DataType> const &
    last_node() { return last_node_; }
    
    inline
    int
    level() { return level_; }
    
    
    
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
    
    ~GraphDynamic() {
        for(int i = 0; i < this->nodes_.size(); i++){
            this->nodes_[i]->unset_connections();
        }
    }
   
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
    
    ~GraphStatic() {
        for(int i = 0; i < this->nodes_.size(); i++){
            this->nodes_[i]->unset_connections();
        }
    }
    
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
            parent_pos = check_pos_is_valid(parent, parent_pos);
            child_pos = check_pos_is_valid(n, child_pos);
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
    connect(
        int i,
        int j,
        int i_pos,
        int j_pos) {
        
        auto n1 = this->get_node(i);
        auto n2 = this->get_node(j);
        i_pos = check_pos_is_valid(n1, i_pos);
        j_pos = check_pos_is_valid(n2, j_pos);
        auto c = std::make_shared<GraphConnection<DataType>>(n1, n2, i_pos, j_pos);
        n1->add_connection(c, i_pos);
        n2->add_connection(c, j_pos);
        this->connections_.push_back(c);
    }
    
    inline
    int
    check_pos_is_valid(
        GraphNodeOP<DataType> const & n,
        int & pos) {
        
        if(pos == -1) {
            Ints avail_pos = n->available_children_pos();
            if(avail_pos.size() == 0) {
                throw GraphException("cannot add connection to node, has not available ends");
            }
            return avail_pos[0];
        }
        
        else {
            if(n->available_pos(pos) == 0) {
                throw GraphException("graph pos is not availabe");
            }
            return pos;
        }
        
    }
    
    inline
    Ints
    get_available_pos(
        GraphNodeOP<DataType> const & n,
        int & pos) {
        
        if(pos == -1) {
            return n->available_children_pos();
        }
        else {
            if(n->available_pos(pos) == 0) {
                std::stringstream ss;
                ss << "graph pos is not available " << pos << std::endl;
                throw GraphException(ss.str());
            }
            Ints r(1);
            r[0] = pos;
            return r;
        }
        
    }
    
    inline
    void
    remove_node(
        int pos) {
        
        auto n = this->get_node(pos);
        for(auto c : n->connections()) {
            if(c == nullptr) { continue; }
            c->disconnect();
            this->connections_.erase(std::remove(this->connections_.begin(), this->connections_.end(),
                                                 c), this->connections_.end());
            
        }
        
        this->nodes_.erase(std::remove(this->nodes_.begin(), this->nodes_.end(),
                                             n), this->nodes_.end());
        
        if(this->nodes_.size() != 0) {
            this->last_node_ = this->nodes_.back();
        }
        else { this->last_node_ = nullptr; }
    }
    
};

#endif /* defined(__RNAMake__graph__) */
