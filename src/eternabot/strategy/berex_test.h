//
//  berex_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__berex_test__
#define __RNAMake__berex_test__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class BerexTest: public Strategy {
public:
    BerexTest() {
        params_ = std::vector<float>(12);
        params_[0] = 0.233331320365;
        params_[1] = 67.1928215814;
        params_[2] = 0.0846125555161;
        params_[3] = 97.8926617326;
        params_[4] = 0.124731583984;
        params_[5] = 76.3152913746;
        params_[6] = -64.8000139266;
        params_[7] = -25.1005947191;
        params_[8] = 1.2097925221;
        params_[9] = 33.3729652592;
        params_[10] = 171.705337159;
        params_[11] = 1.28554296532;
        mean_ = 84.0125821249;
        stdev_ = 8.91633847502;
    }
    
    ~BerexTest() {}
    
public:
    float
    score(FeaturesOP const & features) {
      
        float score = 100;
        score -= fabsf(features->g_count / features->length - params_[0]) * params_[1];
        score -= fabsf(features->u_count / features->length - params_[2]) * params_[3];
        score -= fabsf(features->c_count / features->length - params_[4]) * params_[5];
        
        if     (features->fe < params_[6]) {
            score -= fabsf(features->fe - params_[6]) * params_[8];
        }
        else if(features->fe > params_[7]) {
            score -= fabsf(features->fe - params_[7]) * params_[8];
        }
        
        if     (features->meltpoint < params_[9]) {
            score -= fabsf(features->meltpoint - params_[9]) * params_[11];
        }
        else if(features->meltpoint > params_[10]) {
            score -= fabsf(features->meltpoint - params_[10]) * params_[11];
        }
        
        return score;
    }

};

}


#endif /* defined(__RNAMake__berex_test__) */
