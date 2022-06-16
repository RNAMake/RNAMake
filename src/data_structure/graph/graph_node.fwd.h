//
//  graph_node.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 6/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_graph_node_fwd_h
#define RNAMake_graph_node_fwd_h

#include <memory>
#include <vector>

namespace data_structure {
namespace graph {

enum class GraphNodeType { GraphNodeTypeStatic, GraphNodeTypeDynamic };

template <typename DataType> class GraphNode;

template <typename DataType> class GraphConnection;

template <typename DataType>
using GraphNodeOP = std::shared_ptr<GraphNode<DataType>>;
template <typename DataType>
using GraphConnectionOP = std::shared_ptr<GraphConnection<DataType>>;
template <typename DataType>
using GraphNodeOPs = std::vector<GraphNodeOP<DataType>>;
template <typename DataType>
using GraphConnectionOPs = std::vector<GraphConnectionOP<DataType>>;

} // namespace graph
} // namespace data_structure

#endif
