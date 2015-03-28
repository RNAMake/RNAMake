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
    
    ss_and_seq_ = sstree.get_ss_and_seq();
    if (cofold_) { fr_ = v_.vcofold(ss_and_seq_->seq); }
    else         { fr_ = v_.vfold(ss_and_seq_->seq);   }

    //repopulate data for strategies
    design_data_.fe = fr_.free_energy;
    design_data_.gc = sstree.gc_pairs();
    design_data_.gu = sstree.gu_pairs();
    design_data_.ua = sstree.ua_pairs();
    design_data_.dotplot = v_.dotplot(ss_and_seq_->seq);
    if(design_data_.pairmap.size() == 0) {
        design_data_.pairmap = sstree.get_pairmap();
        design_data_.length = (int)ss_and_seq_->seq.length();
    }
    design_data_.sstree = sstree;
    design_data_.seq = ss_and_seq_->seq;
    
    float total_score = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        total_score += ((s->score(design_data_) - s->mean()) / (s->stdev())) * weights_[i];
        i++;
    }
        
    return  (total_score * stdev_) + mean_ ;
    
}