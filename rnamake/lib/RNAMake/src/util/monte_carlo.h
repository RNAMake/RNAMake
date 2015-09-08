//
//  monte_carlo.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__monte_carlo__
#define __RNAMake__monte_carlo__

#include <stdio.h>
#include <math.h>

//RNAMake Headers
#include "util/random_number_generator.h"

class MonteCarlo {
public:
    inline
    MonteCarlo(
        float temperature=1.0f):
    temperature_(temperature),
    rng_(RandomNumberGenerator()){}
    
    ~MonteCarlo() {}

    inline
    int
    accept(
        float current,
        float next) {
        
        if(next < current) { return 1; }
        
        score_ = exp((current - next) / temperature_);
        if(rng_.rand() < score_) { return 1; }
        
        return 0;
    }
    
private:
    float temperature_;
    float score_;
    RandomNumberGenerator rng_;
    
};

#endif /* defined(__RNAMake__monte_carlo__) */
