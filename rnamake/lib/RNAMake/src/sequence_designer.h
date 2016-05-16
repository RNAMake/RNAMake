//
//  sequence_designer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 2/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sequence_designer__
#define __RNAMake__sequence_designer__

#include <stdio.h>
#include "FileIO.h"
#include "random_number_generator.h"
#include "secondary_structure_tree.h"
#include "sequence_design_scorer.h"
#include "eternabot_scorer.h"

struct SeqandScore {
    String seq;
    float score;
};

struct seq_and_score_less_than_key {
    inline bool operator() (
        SeqandScore const & ss1,
        SeqandScore const & ss2) {
        
        return (ss1.score > ss2.score);
    }
};

typedef std::vector<SeqandScore> SeqandScores;

class SequenceDesigner {
public:
    SequenceDesigner(){
        rng_  = RandomNumberGenerator();
        basepairs_ = split_str_by_delimiter("GC,CG,AU,UA", ",");
        steps_ = 100;
        scorer_ = SequenceDesignScorerOP( new EternabotScorer() );
        score_ = 10000;
        results_ = SeqandScores(steps_/10);
        all_results_ = SeqandScores(steps_);
    }
    
    ~SequenceDesigner() {}
    
public:
    
    String
    design(
        String const &,
        String const &);
    
    void
    mutate_basepair(
        SecondaryStructureNodeOP & bp_node);
    
    
public: //getters
    
    inline
    float
    score() { return score_; }
    
    inline
    std::vector<float> const &
    scores() { return scorer_->scores(); }
    
    inline
    SeqandScores const &
    top_seqs() { return results_; }
    
    
private:
    RandomNumberGenerator rng_;
    SecondaryStructureTree ss_tree_;
    SequenceDesignScorerOP scorer_;
    Strings basepairs_;
    String last_bp_type_, cseq_;
    SeqandScores results_, all_results_;
    float score_, best_score_;
    int steps_;
};

#endif /* defined(__RNAMake__sequence_designer__) */