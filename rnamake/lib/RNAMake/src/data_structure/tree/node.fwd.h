//
//  node.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_node_fwd_h
#define RNAMake_node_fwd_h

#include <memory>
#include <vector>

template <class DataType>
class Node;

template <class DataType>
using NodeOP  = std::shared_ptr<Node<DataType>>;

template <class DataType>
using NodeOPs = std::vector<NodeOP<DataType>>;


#endif
