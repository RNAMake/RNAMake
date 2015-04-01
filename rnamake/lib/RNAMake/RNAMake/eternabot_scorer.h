//
//  eternabot_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 3/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__eternabot_scorer__
#define __RNAMake__eternabot_scorer__

#include <stdio.h>
#include <map>
#include "secondary_structure_node.h"
#include "secondary_structure_tree.h"
#include "sequence_design_scorer.h"
#include "eterna_bot_strategy.h"


class EternabotScorer : public SequenceDesignScorer {
public:
    EternabotScorer() :
    design_data_ ( EternabotScorerData() ),
    strategies_ ( EternabotStrategyOPs() ),
    weights_ ( std::vector<float>() ) {
        strategies_.push_back( EternabotStrategyOP ( new ABasicTest() ) );
        strategies_.push_back( EternabotStrategyOP ( new ClearPlotStackCapsandSafeGC() ) );
        strategies_.push_back( EternabotStrategyOP ( new BerexTest() ) );
        strategies_.push_back( EternabotStrategyOP ( new NumofYellowNucleotidesperLengthofString() ) );
        strategies_.push_back( EternabotStrategyOP ( new DirectionofGCPairsinMultiLoops() ) );
        
        weights_.push_back(-0.0981965481737);
        weights_.push_back(0.44654228578);
        weights_.push_back(0.325776127974);
        weights_.push_back(0.173103278813);
        weights_.push_back(0.139425514791);
        
        mean_ = 84.8005952381;
        stdev_ = 16.4725276237;

    }
    
    ~EternabotScorer() {}
    
public:

    float
    score_sstree(
        SecondaryStructureTree const &);

private:
    EternabotScorerData design_data_;
    EternabotStrategyOPs strategies_;
    std::vector<float> weights_;
    float mean_, stdev_;
    
};

#endif /* defined(__RNAMake__eternabot_scorer__) */




