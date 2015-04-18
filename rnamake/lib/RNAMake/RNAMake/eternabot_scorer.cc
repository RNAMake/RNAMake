//
//  eternabot_scorer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 3/21/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "eternabot_scorer.h"

float
EternabotScorer::score_sstree(
    SecondaryStructureTree const & sstree) {

    String seq = sstree.seq();
    
    //repopulate data for strategies
    design_data_.gc = sstree.gc_pairs();
    design_data_.gu = sstree.gu_pairs();
    design_data_.ua = sstree.ua_pairs();
    design_data_.dotplot = v_.bp_probabilities(seq);
    design_data_.fe = v_.free_energy();
    
    if(design_data_.pairmap.size() == 0) {
        design_data_.pairmap = sstree.get_pairmap();
        design_data_.length = (int)seq.length();
    }
    
    design_data_.sstree = sstree;
    design_data_.seq = seq;
        
    float total_score = 0;
    float score = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        score = s->score(design_data_);
        scores_[i] = score;
        total_score += ((score - s->mean()) / (s->stdev())) * weights_[i];
        i++;
    }
        
    return  (total_score * stdev_) + mean_ ;
    
}