//
//  sequence_designer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "sequence_designer.h"

String
SequenceDesigner::design(
    String const & sequence,
    String const & structure) {
    
    scorer_->setup(structure);
    ss_tree_ = SecondaryStructureTree(structure, sequence);
    SecondaryStructureNodeOPs designable_bps = ss_tree_.get_designable_bps();
    for(auto & bp_node : designable_bps) { mutate_basepair(bp_node); }
    
    int bps_pos = 0, bps_size = (int)designable_bps.size();
    float new_score;
    score_ = scorer_->score_sstree(ss_tree_);
    while(isnan(score_)) {
        for(auto & bp_node : designable_bps) { mutate_basepair(bp_node); }
        score_ = scorer_->score_sstree(ss_tree_);
    }
    ss_and_seq_ = ss_tree_.get_ss_and_seq();
    cseq_ = ss_and_seq_->seq;
    best_score_ = -1000;
    steps_ = 10;
    float T = 2;
    float diceroll;
    float prob;
    int i = 0;
    int count = 0;
    while(i < steps_) {
        bps_pos = rng_.randrange(bps_size);
        std::cout << bps_size << " " << bps_pos << std::endl;
        mutate_basepair(designable_bps[bps_pos]);
        new_score = scorer_->score_sstree(ss_tree_);
        if(isnan(new_score) || new_score < -100000) {
            count++;
            if(count > 100) {
                break;
            }
            designable_bps[bps_pos]->bp_type(last_bp_type_);
            continue;
        }
        count = 0;
        i++;
        if(new_score > score_) {
            score_ = new_score;
            if(best_score_ < score_) {
                ss_and_seq_ = ss_tree_.get_ss_and_seq();
                cseq_ = ss_and_seq_->seq;
                best_score_ = score_;
            }
            continue;
        }
        
        prob = expf((new_score - score_) / T);
        diceroll = rng_.rand();
        if(diceroll < prob) {
            score_ = new_score;
            continue;
        }
        
        if(i % 100) { T = T*0.8; }
        
        designable_bps[bps_pos]->bp_type(last_bp_type_);
    }
    
    score_ = best_score_;
    
    return cseq_;
    
}


void
SequenceDesigner::mutate_basepair(
    SecondaryStructureNodeOP & bp_node) {
    
    int pos = rng_.randrange(4);
    last_bp_type_ = bp_node->bp_type();
    bp_node->bp_type(basepairs_[pos]);
    
}






