//
//  scorer.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot/scorer.h"

namespace eternabot {


void
Scorer::setup(secondary_structure::PoseOP const & p) {
    features_ = generator_.get_features(p);
    
}

float
Scorer::score_secondary_structure(secondary_structure::PoseOP const & p) {
    generator_.update_features(features_, p);
    
    total_score_ = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        scores_[i] = s->score(features_);
        total_score_ += scores_[i]*weights_[i];
        //std::cout << i << " " << scores_[i] << " " << weights_[i] << std::endl;
        i++;
    }
    //std::cout << total_score_ << std::endl;
    //exit(0);
    return total_score_;
    //return (total_score_ * stdev_) + mean_ ;
}

float
Scorer::print_scores(secondary_structure::PoseOP const & p) {
    generator_.update_features(features_, p);

    total_score_ = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        scores_[i] = s->score(features_);
        total_score_ += scores_[i]*weights_[i];
        std::cout << s->name() << " " << scores_[i] << " " << weights_[i] << std::endl;
        i++;
    }
    return total_score_;
}

Floats
Scorer::get_scores(secondary_structure::PoseOP const & p) {
    generator_.update_features(features_, p);
    auto scores = Floats();

    total_score_ = 0;
    int i = 0;
    for(auto const & s : strategies_) {
        scores.push_back(s->score(features_));


    }
    return scores;
}


}