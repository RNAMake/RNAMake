//
//  scorer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "scorer.h"

namespace eternabot {


void
Scorer::setup(sstruct::PoseOP const & p) {
    features_ = generator_.get_features(p);
    
}

float
Scorer::score_secondary_structure(sstruct::PoseOP const & p) {
    generator_.update_features(features_, p);
    
    total_score_ = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        scores_[i] = s->score(features_);
        //total_score_ += ((scores_[i] - s->mean()) / (s->stdev())) * weights_[i];
        total_score_ += scores_[i]*weights_[i];
        //std::cout << scores_[i] << " " << weights_[i] << " " << scores_[i]*weights_[i] << " " <<  ((scores_[i] - s->mean()) / (s->stdev())) * weights_[i] << std::endl;
        i++;
    }
    //std::cout << total_score_ << std::endl;
    //exit(0);
    return total_score_;
    //return (total_score_ * stdev_) + mean_ ;
}


}