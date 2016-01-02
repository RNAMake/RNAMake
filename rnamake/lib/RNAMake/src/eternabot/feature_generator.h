//
//  feature_generator.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/1/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__feature_generator__
#define __RNAMake__feature_generator__

#include <stdio.h>

#include "vienna/vienna.h"
#include "secondary_structure/pose.h"

namespace eternabot {

struct Features {
public:
    
    Features():
    gu(0), gc(0), ua(0),
    meltpoint(97), fe(0),
    pairmap ( std::map<int, int>() )
    {}
    
    
    float gu, gc, ua;
    float meltpoint, fe;
    vienna::plists dotplot;
    std::map<int, int> pairmap;
    
};

typedef std::shared_ptr<Features> FeaturesOP;
    
class FeatureGenerator {
public:
    FeatureGenerator():
    v_(vienna::Vienna()) {}
    
    ~FeatureGenerator() {}
    
public:
    
    FeaturesOP
    get_features(
        sstruct::PoseOP const &);
    
    void
    update_features(
        FeaturesOP &,
        sstruct::PoseOP const &);
  
private:
    vienna::Vienna v_;
    
};
    
}

#endif /* defined(__RNAMake__feature_generator__) */
