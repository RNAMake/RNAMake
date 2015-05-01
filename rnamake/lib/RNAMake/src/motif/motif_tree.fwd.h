//
//  motif_tree.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/4/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_motif_tree_fwd_h
#define RNAMake_motif_tree_fwd_h
#include <vector>
#include <memory>

class MotifTreeConnection;
typedef std::vector<MotifTreeConnection> MotifTreeConnections;
typedef std::shared_ptr<MotifTreeConnection> MotifTreeConnectionOP;
typedef std::vector<MotifTreeConnectionOP> MotifTreeConnectionOPs;
class MotifTreeNode;
typedef std::shared_ptr<MotifTreeNode> MotifTreeNodeOP;
typedef std::vector<MotifTreeNodeOP> MotifTreeNodeOPs;
class MotifTree;


#endif
