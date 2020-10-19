//
//  berex_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake_mod_berex_test__
#define __RNAMake_mod_berex_test__

#include <stdio.h>

#include "eternabot/strategy.h"

namespace eternabot {

class ModifiedBerexTest: public Strategy {
public:
    ModifiedBerexTest() {
        params_ = std::vector<float>(12);
        params_[0] = 0.150150639777;
        params_[1] = 116.231035752;
        params_[2] = 0.068678024689;
        params_[3] = 129.231308029;
        params_[4] = 0.195823646983;
        params_[5] = 117.625896452;
        params_[6] = -68.3488778805;
        params_[7] = -30.2069356813;
        params_[8] = 1.12021566776;
        params_[9] = 35.4749400799;
        params_[10] = 133.662333848;
        params_[11] = 1.36225265987;
        mean_ = 84.0125821249;
        stdev_ = 8.91633847502;
        name_ = "ModifiedBerexTest";
    }
    
    ~ModifiedBerexTest() {}
    
public:
    float
    score(FeaturesOP const & features) {
      
        float score = 100;

        if(features->length > 30) {
            score -= abs(float(features->g_count) / float(features->length) - params_[0]) * params_[1];
            score -= abs(float(features->u_count) / float(features->length) - params_[2]) * params_[3];
            score -= abs(float(features->c_count) / float(features->length) - params_[4]) * params_[5];
        }

        float weight = exp(-((features->length - 100)*(features->length - 100))/5000)/(sqrt(2*3.14*1))*2.5;
        if     (features->fe < params_[6]) {
           score -= abs(features->fe - params_[6]) * params_[8] * weight;
        }
        else if(features->fe > params_[7]) {
           score -= abs(features->fe - params_[7]) * params_[8] * weight;
        }

        if     (features->meltpoint < params_[9]) {
            score -= abs(features->meltpoint - params_[9]) * params_[11] * weight;
        }
        else if(features->meltpoint > params_[10]) {
            score -= abs(features->meltpoint - params_[10]) * params_[11] * weight;
        }

        /*if     (features->fe < params_[6]) {
            score -= abs(features->fe - params_[6]) * params_[8];
        }
        else if(features->fe > params_[7]) {
            score -= abs(features->fe - params_[7]) * params_[8];
        }*/
        
        /*if     (features->meltpoint < params_[9]) {
            score -= abs(features->meltpoint - params_[9]) * params_[11];
        }
        else if(features->meltpoint > params_[10]) {
            score -= abs(features->meltpoint - params_[10]) * params_[11];
        }*/
        
        return score;
    }

};

}


#endif /* defined(__RNAMake__berex_test__) */
