//
//  mtst_to_met_converter.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/20/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__mtst_to_met_converter__
#define __RNAMake__mtst_to_met_converter__

#include "motif_tree_state_tree.h"
#include "motif_ensemble_tree.fwd.h"
#include "motif_tree.h"
#include "types.h"
#include "pose.h"
#include <stdio.h>
#include <map>

class MTSTtoMETConverter {
public:
    MTSTtoMETConverter() {
        node_num_map_ = std::map<int, int>();
    }
    
    ~MTSTtoMETConverter() {}

public:
    MotifEnsembleTreeOP
    convert(
        MotifTreeStateTree const &,
        int start_pos = 1);
    
private:
    
    ChainOP
    _get_start_chain(
        MotifTreeNodeOP const &);
    
    String
    _get_next_bp(
        MotifTreeNodeOP const &,
        MotifTreeNodeOP const &,
        ChainOP const &,
        int flip = 0);
    
    String
    _add_helix(
        MotifTreeNodeOP const &,
        ChainOP const &,
        String const &,
        int const &,
        int const &,
        MotifEnsembleTreeNodeOP const &,
        int const &);
    
    int
    _get_chain_pos(
        ResidueOP const &);
    
private:
    PoseOP p_;
    MotifTree mt_;
    MotifEnsembleTreeOP met_;
    String dseq_;
    std::map<int, int> node_num_map_;
    
    
};

#endif /* defined(__RNAMake__mtst_to_met_converter__) */
