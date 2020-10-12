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

namespace util {

class MonteCarlo {
public:
    inline
    MonteCarlo(
            float temperature = 1.0f) :
            temperature_(temperature),
            rng_(RandomNumberGenerator()) {}

    MonteCarlo(
            float temperature,
            int seed
            ) :
            temperature_(temperature),
            rng_(RandomNumberGenerator(seed)) {}

    ~MonteCarlo() {}

    //inline
    //int
    //accept(
    //        float current,
    //        float next) {

    //    if (next < current) { return 1; }

    //    score_ = exp((current - next) / temperature_);
    //    if (rng_.rand() < score_) { return 1; }

    //    return 0;
    //}
    template<class T, class U, typename std::enable_if<std::is_floating_point<T>::value&std::is_floating_point<U>::value>::type* = nullptr>
    inline
    int
    accept(
            T current,
            U next) {

        if (next < current) { return 1; }

        score_ = exp((current - next) / temperature_);
        if (rng_.rand() < score_) { return 1; }

        return 0;
    }
public:
    inline
    void
    set_temperature(
            float new_temp) { temperature_ = new_temp; }

    inline
    float
    get_temperature() { return temperature_; }


    inline
    void
    scale_temperature(
            float scale) {
        temperature_ *= scale;
    }

private:
    float temperature_;
    float score_;
    RandomNumberGenerator rng_;

};

}

#endif /* defined(__RNAMake__monte_carlo__) */
