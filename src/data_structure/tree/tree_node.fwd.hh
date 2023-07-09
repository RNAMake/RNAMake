//
//  tree_node.fwd.hh
//  RNAMake
//
//  Created by Joseph Yesselman on 7/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_tree_node_fwd_hh
#define RNAMake_tree_node_fwd_hh

namespace data_structure {
namespace tree {

template<class DataType>
class TreeNode;

template<class DataType>
using TreeNodeOP  = std::shared_ptr<TreeNode<DataType>>;

template<class DataType>
using TreeNodeOPs = std::vector<TreeNodeOP<DataType>>;

}
}

#endif
