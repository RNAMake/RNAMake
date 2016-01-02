//
//  feature_generator.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 1/1/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#include "eternabot/feature_generator.h"


namespace eternabot {

FeaturesOP
FeatureGenerator::get_features(
    sstruct::PoseOP const & p) {
    
    auto features = std::make_shared<Features>();
    for(auto const & bp : p->basepairs()) {
        features->pairmap[bp->res1()->num()] = bp->res2()->num();
        features->pairmap[bp->res2()->num()] = bp->res1()->num();
    }
    
    update_features(features, p);
    return features;
    
}

void
FeatureGenerator::update_features(
    FeaturesOP & features,
    sstruct::PoseOP const & p) {
    
    for(auto const & bp : p->basepairs()) {
        if     (bp->res1()->res_type() == 1 && bp->res2()->res_type() == 2) { features->gc += 1; }
        else if(bp->res1()->res_type() == 2 && bp->res2()->res_type() == 1) { features->gc += 1; }
        else if(bp->res1()->res_type() == 0 && bp->res2()->res_type() == 3) { features->ua += 1; }
        else if(bp->res1()->res_type() == 3 && bp->res2()->res_type() == 0) { features->ua += 1; }
        else if(bp->res1()->res_type() == 2 && bp->res2()->res_type() == 3) { features->gu += 1; }
        else if(bp->res1()->res_type() == 3 && bp->res2()->res_type() == 2) { features->gu += 1; }
    }
    
    features->dotplot = v_.bp_probabilities(p->sequence());
    features->fe = v_.free_energy();
    
    
}

}