//
//  sequence_optimizer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/11/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_optimizer__
#define __RNAMake__sequence_optimizer__

#include <stdio.h>
#include "eternabot/sequence_designer.h"
#include "motif_data_structures/motif_graph.h"
#include "motif_data_structures/motif_tree.h"

struct SequenceOptimizerResult {
    inline
    SequenceOptimizerResult(
        MotifTreeOP const & nmt,
        float nscore):
    motif_tree(nmt),
    score(nscore)
    {}
    
    MotifTreeOP motif_tree;
    float score;
};

typedef std::shared_ptr<SequenceOptimizerResult> SequenceOptimizerResultOP;

class SequenceOptimizer {
public:
    SequenceOptimizer() {}
    
    ~SequenceOptimizer() {}
    
public:
    
    SequenceOptimizerResultOP
    optimize(
        MotifGraphOP &,
        int,
        int,
        int,
        int);
    
private:
    MotifTreeOP mt_;
    eternabot::SequenceDesigner designer_;
    eternabot::SequenceDesignerResultOPs designer_results_;
};

#endif /* defined(__RNAMake__sequence_optimizer__) */
