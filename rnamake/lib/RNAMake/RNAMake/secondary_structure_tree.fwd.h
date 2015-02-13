//
//  secondary_structure_tree.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/11/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_secondary_structure_tree_fwd_h
#define RNAMake_secondary_structure_tree_fwd_h
#include <vector>

class SecondaryStructureNode;
typedef std::shared_ptr<SecondaryStructureNode> SecondaryStructureNodeOP;
typedef std::vector<SecondaryStructureNodeOP> SecondaryStructureNodeOPs;
typedef std::vector<SecondaryStructureNode*> SecondaryStructureNodes;

#endif
