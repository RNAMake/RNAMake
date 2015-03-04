//
//  sequence_design_scorer.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 3/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "sequence_design_scorer.h"

float
SequenceDesignScorer::score(
    String const & seq) {
    
    if (cofold_) { fr_ = v_.vcofold(seq); }
    else         { fr_ = v_.vfold(seq);   }
        
    int diff = 0;
    for (int i = 0; i < ss_.length(); i++) {
        if(ss_[i] != fr_.structure[i]) { diff++; }
    }
    
    score_ = diff*10 + fr_.free_energy;
    return score_;
    
}