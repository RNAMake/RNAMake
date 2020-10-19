//
//  a_basic_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__mod_basic_test__
#define __RNAMake__mod_basic_test__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class ModifiedABasicTest : public Strategy {
public:
    ModifiedABasicTest() {
        params_ = Floats(7);
        params_[0] = 0.303497971269;
        params_[1] = 92.9893755247;
        params_[2] = -1.37878787864;
        params_[3] = 0.512804062262;
        params_[4] = 0.477932936507;
        params_[5] = 84.4793979751;
        params_[6] = 124.345433009;
        mean_ = 83.5007560083;
        stdev_ = 10.5290224709;
        
    }
    
    ~ModifiedABasicTest() {}
    
    inline
    float
    score(FeaturesOP const & features) {
        float total_pairs = features->gc + features->gu + features->ua;
        float score = 100;
        if(total_pairs > 0) {
            score -= fabs(features->ua / total_pairs - params_[0]) * params_[1];
        }
        float target_fe = params_[2] * (total_pairs+1);
        score -= fabs(target_fe - features->fe) * params_[3];

        if(features->meltpoint < params_[5]) {
            score -= fabs(features->meltpoint - params_[5]) * params_[4];
        }
        else if(features->meltpoint > params_[6]) {
            score -= fabs(features->meltpoint - params_[6]) * params_[4];
        }
        return score;
    }
    
};
    
}

#endif /* defined(__RNAMake__a_basic_test__) */
