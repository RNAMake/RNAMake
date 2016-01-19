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
#include <algorithm>
#include <queue>
#include <map>
#include <cassert>

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
    
    virtual
    ~Graph() {
        for(int i = 0; i < nodes_.size(); i++){
            nodes_[i]->unset_connections();
        }
    }
    
public:
    typedef GraphIterator<DataType> iterator;
    typedef const GraphIterator<DataType> const_iterator;
    friend class GraphIterator<DataType>;
    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;
    
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
    GraphNodeOP<DataType>
    oldest_node() {
        auto node = last_node_;
        assert(node != nullptr && "attemped to call oldest_node but there are no nodes");
        
        for(auto const & n : nodes_) {
            if(n->index() < node->index()) { node = n; }
        }
        
        return node;
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
    nodes() const {
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
    GraphNodeOPs<DataType> empty_;
    int index_, level_;
};

template <typename DataType>
class GraphIterator {
    friend class Graph<DataType>;
    
public:
    GraphIterator() {}
    
public:
    GraphIterator & operator++ ();
    GraphIterator operator++ (int);
    
    bool operator== (const GraphIterator& rhs) const;
    bool operator!= (const GraphIterator& rhs) const;
    GraphNodeOP<DataType> & operator* ();

    
private:
    GraphNodeOPs<DataType> nodes_, leafs_;
    GraphNodeOP<DataType> current_;
    std::queue<GraphNodeOP<DataType>> queue_;
    std::map<GraphNodeOP<DataType>, int> seen_;
    
    GraphIterator(
        GraphNodeOPs<DataType> const & nodes):
    nodes_(nodes)
    {
        if(nodes_.size() == 0) {
            current_ = nullptr;
            return;
        }
        
        queue_ = std::queue<GraphNodeOP<DataType>>();
        seen_  = std::map<GraphNodeOP<DataType>, int> ();
        
        int active_conn = 0;
        for(auto const & n : nodes) {
            active_conn = 0;
            for (auto const & c : n->connections()) {
                if(c != nullptr) { active_conn += 1; }
            }
            
            if(active_conn < 2) {
                leafs_.push_back(n);
            }
        }
                
        if(leafs_.size() > 0) {
            current_ = leafs_[0];
        }
        else {
            current_ = nodes_[0];
        }
        
        seen_[current_] = 1;
        
    }
};




template <typename DataType>
typename Graph<DataType>::iterator
Graph<DataType>::begin()  {
    return iterator(nodes_);
}

template <typename DataType>
typename Graph<DataType>::iterator
Graph<DataType>::end()  {
    return iterator(empty_);
}

template <typename DataType>
typename Graph<DataType>::const_iterator
Graph<DataType>::begin() const {
    return iterator(nodes_);
}

template <typename DataType>
typename Graph<DataType>::const_iterator
Graph<DataType>::end()  const {
    return iterator(empty_);
}

template <typename DataType>
GraphNodeOP<DataType>&
GraphIterator<DataType>::operator* () {
    return current_;
}

template <typename DataType>
GraphIterator<DataType>&
GraphIterator<DataType>::operator++() {
    for(auto const & c : current_->connections()) {
        if(c == nullptr) { continue; }
        auto n = c->partner(current_->index());
        if(seen_.find(n) != seen_.end()) {
            continue;
        }
        queue_.push(n);
        seen_[n] = 1;
        
    }
    
    if(!queue_.empty()) {
        current_ = queue_.front();
        queue_.pop();
    }
    else {
        int found = 0;
        for(auto const & n : leafs_) {
            if(seen_.find(n) == seen_.end()) {
                seen_[n] = 1;
                current_ = n;
                found = 1;
                break;
            }
        }
        
        if(!found) {
            current_ = nullptr;
        }
    }
    
    
    return *this;
}

template <typename DataType>
bool
GraphIterator<DataType>::operator== (
    GraphIterator<DataType> const & rhs) const {
    return current_ == rhs.current_;
}

template <typename DataType>
bool
GraphIterator<DataType>::operator!= (
    GraphIterator<DataType> const& rhs) const {
    return current_ != rhs.current_;
}




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
    
    inline
    void
    connect(
        int i,
        int j) {
        
        auto n1 = this->get_node(i);
        auto n2 = this->get_node(j);
        auto c = std::make_shared<GraphConnection<DataType>>(n1, n2, 0, 0);
        n1->add_connection(c);
        if(n1 != n2) { n2->add_connection(c); }
        this->connections_.push_back(c);
        
        
    }
    
    
};

template <typename DataType>
class GraphStatic : public Graph<DataType> {
public:
    inline
    GraphStatic(): Graph<DataType>() {}
    
    inline
    GraphStatic(
        GraphStatic<DataType> const & g) {
        this->nodes_ = GraphNodeOPs<DataType>(g.nodes_.size());
        int i = 0;
        for(auto const & n : g.nodes_) {
            this->nodes_[i] = std::make_shared<GraphNodeStatic<DataType>>(*n);
            i++;
        }
        
        i = 0;
        int j, ei, ej;
        for(auto const & c : g.connections_) {
            i = c->node_1()->index();
            j = c->node_2()->index();
            ei = c->end_index(i);
            ej = c->end_index(j);
            connect(i, j, ei, ej);
        }
        
        if(g.last_node_ != nullptr) {
            this->last_node_ = this->get_node(g.last_node_->index());
        }
        
        this->level_ = g.level_;
        this->index_ = g.index_;
        
    }
    
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
        int n_children = 0,
        int orphan = 0) {
        
        GraphNodeOP<DataType> parent = this->last_node_;
        auto n = std::make_shared<GraphNodeStatic<DataType>>(data, this->index_, this->level_,
                                                             n_children);

        if(parent_index != -1) { parent = this->get_node(parent_index); }
        if(orphan == 1) {
            parent = nullptr;
        }
        if(parent != nullptr) {
            parent_pos = check_pos_is_valid(parent, parent_pos);
            child_pos = check_pos_is_valid(n, child_pos);
            auto c = std::make_shared<GraphConnection<DataType>>(parent, n, parent_pos, child_pos);
            parent->add_connection(c, parent_pos);
            n->add_connection(c, child_pos);
            this->connections_.push_back(c);
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
            auto partner = c->partner(n->index());
            n->remove_connection(c);
            partner->remove_connection(c);
            this->connections_.erase(std::remove(this->connections_.begin(),
                                                 this->connections_.end(), c));
            c->disconnect();
            
        }
        
        this->nodes_.erase(std::remove(this->nodes_.begin(), this->nodes_.end(),
                                             n), this->nodes_.end());
        
        if(this->nodes_.size() != 0) {
            this->last_node_ = this->nodes_.back();
        }
        else { this->last_node_ = nullptr; }
    }
    
    inline
    void
    remove_level(
        int level) {
        
        int pos = 0;
        while(pos < this->nodes_.size()) {
            auto n = this->nodes_[pos];
            if(n->level() >= level) { remove_node(n->index()); continue; }
            pos++;
        }
        
    }
    
};




























#endif /* defined(__RNAMake__graph__) */
