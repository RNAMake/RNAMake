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
    level_(0),
    node_type_( GraphNodeType::GraphNodeTypeDynamic ){}
    
    Graph(
        GraphNodeType const & type):
    node_type_(type),
    index_(0),
    level_(0)
    {}
    
    ~Graph() {}
    
public:
    
    inline
    int
    add_data(
        DataType const & data,
        int n_children=0,
        int parent_index=-1,
        int parent_pos=-1,
        int child_pos=-1) {
        auto n = std::make_shared<GraphNode<DataType>>(data, index_, level_, node_type_, n_children);
        GraphNodeOP<DataType> parent = last_node_;
        if(parent_index != -1) { parent = get_node(parent_index);}
        if(parent == nullptr) {
            nodes_.push_back(n);
            index_++;
            last_node_ = n;
            return 0;
        }
        else if (parent_pos == -1 && node_type_ == GraphNodeType::GraphNodeTypeStatic) {
            Ints avail_child_pos = parent->available_children_pos();
            if(avail_child_pos.size() == 0) {
                throw GraphException("cannot use node as parent as it has not spots for children\n");
            }
            parent_pos = avail_child_pos[0];
        }
        
        auto connection = std::make_shared<GraphConnection<DataType>>(parent, n, parent_pos, child_pos);
        parent->add_connection(connection, parent_pos);
        n->add_connection(connection, child_pos);
        
        nodes_.push_back(n);
        connections_.push_back(connection);
        index_++;
        last_node_ = n;
        
        return n->index();
        
    }
    
    inline
    GraphNodeOP<DataType> const &
    get_node(
        int index) {
        
        for(auto const & n : nodes_) {
            if(n->index() == index) { return n; }
        }
        
        throw GraphException("cannot find node with index");
    }
    
     

    
private:
    GraphNodeOP<DataType> last_node_;
    GraphNodeOPs<DataType> nodes_;
    GraphConnectionOPs<DataType> connections_;
    int index_, level_;
    GraphNodeType node_type_;
};

template <typename DataType>



#endif /* defined(__RNAMake__graph__) */
