//
//  motif_state_seatch_node.fwd.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMake_motif_state_seatch_node_fwd_h
#define RNAMake_motif_state_seatch_node_fwd_h

#include <memory>
#include <vector>

namespace motif_search {

class MotifStateSearchNode;

typedef std::shared_ptr<MotifStateSearchNode> MotifStateSearchNodeOP;
typedef std::vector<MotifStateSearchNodeOP> MotifStateSearchNodeOPs;

}

#endif
