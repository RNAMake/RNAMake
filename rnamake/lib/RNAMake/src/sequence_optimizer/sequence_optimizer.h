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

struct OptimizedSequence {
    inline
    OptimizedSequence(
        String const & nsequence,
        float const & nclose_distance,
        float const & neternabot_score):
    sequence(nsequence),
    close_distance(nclose_distance),
    eternabot_score(neternabot_score)
    {}
    
    String sequence;
    float close_distance, eternabot_score;
};

typedef std::shared_ptr<SequenceOptimizerResult> SequenceOptimizerResultOP;
typedef std::shared_ptr<OptimizedSequence>       OptimizedSequenceOP;
typedef std::vector<OptimizedSequenceOP>         OptimizedSequenceOPs;


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
    
    OptimizedSequenceOPs
    get_optimized_sequences(
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
