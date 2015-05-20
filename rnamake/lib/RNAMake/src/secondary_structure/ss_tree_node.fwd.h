//
//  ss_tree_node.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_ss_tree_node_fwd_h
#define RNAMake_ss_tree_node_fwd_h
#include <vector>
#include <memory>

class SS_TreeNode;
typedef std::shared_ptr<SS_TreeNode> SS_TreeNodeOP;
typedef std::vector<SS_TreeNodeOP>   SS_TreeNodeOPs;

struct SS_TreeNodeProto;
typedef std::shared_ptr<SS_TreeNodeProto> SS_TreeNodeProtoOP;
typedef std::vector<SS_TreeNodeProtoOP>   SS_TreeNodeProtoOPs;

#endif
