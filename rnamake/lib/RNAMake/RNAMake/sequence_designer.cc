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
        std::cout << score_ << std::endl;
    }
    ss_and_seq_ = ss_tree_.get_ss_and_seq();
    cseq_ = ss_and_seq_->seq;
    steps_ = 1000;
    
    for(int i = 0; i < steps_; i++) {
        bps_pos = rng_.randrange(bps_size);
        mutate_basepair(designable_bps[bps_pos]);
        ss_and_seq_ = ss_tree_.get_ss_and_seq();
        new_score = scorer_->score_sstree(ss_tree_);
        if(new_score > score_) {
            score_ = new_score;
            cseq_ = ss_and_seq_->seq;
        }
        else {
            designable_bps[bps_pos]->bp_type(last_bp_type_);
        }
    }
    
    return cseq_;
    
}


void
SequenceDesigner::mutate_basepair(
    SecondaryStructureNodeOP & bp_node) {
    
    int pos = rng_.randrange(4);
    last_bp_type_ = bp_node->bp_type();
    bp_node->bp_type(basepairs_[pos]);
    
}






