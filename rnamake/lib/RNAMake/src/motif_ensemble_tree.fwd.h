//
//  motif_ensemble_tree.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_motif_ensemble_tree_fwd_h
#define RNAMake_motif_ensemble_tree_fwd_h
#include <vector>
#include <memory>

class MotifEnsembleTreeNode;
typedef std::shared_ptr<MotifEnsembleTreeNode> MotifEnsembleTreeNodeOP;
typedef std::vector<MotifEnsembleTreeNodeOP> MotifEnsembleTreeNodeOPs;

class MotifEnsembleTree;
typedef std::shared_ptr<MotifEnsembleTree> MotifEnsembleTreeOP;


#endif
