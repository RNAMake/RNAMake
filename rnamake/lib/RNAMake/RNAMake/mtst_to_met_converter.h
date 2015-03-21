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
#include "pose.h"
#include <stdio.h>

class MTSTtoMETConverter {
public:
    MTSTtoMETConverter() {}
    
    ~MTSTtoMETConverter() {}

public:
    MotifEnsembleTreeOP
    convert(
        MotifTreeStateTree const &);
    
private:
    
    ChainOP
    _get_start_chain(
        MotifTreeNodeOP const & n);
    
private:
    PoseOP p_;
    MotifTree mt_;
    String dseq_;
    
    
};

#endif /* defined(__RNAMake__mtst_to_met_converter__) */
