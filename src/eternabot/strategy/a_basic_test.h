//
//  a_basic_test.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__a_basic_test__
#define __RNAMake__a_basic_test__

#include <stdio.h>

#include <eternabot/strategy.h>

namespace eternabot {

class ABasicTest : public Strategy {
public:
  ABasicTest() {
    _params = Reals(7);
    _params[0] = 0.303497971269;
    _params[1] = 92.9893755247;
    _params[2] = -1.37878787864;
    _params[3] = 0.512804062262;
    _params[4] = 0.477932936507;
    _params[5] = 84.4793979751;
    _params[6] = 124.345433009;
    _mean = 83.5007560083;
    _stdev = 10.5290224709;
  }

  ~ABasicTest() {}

  inline float score(FeaturesOP const &features) {
    float total_pairs = features->gc + features->gu + features->ua;
    float score = 100;
    if (total_pairs > 0) {
      score -= fabs(features->ua / total_pairs - _params[0]) * _params[1];
    }
    float target_fe = _params[2] * total_pairs;
    score -= fabs(target_fe - features->fe) * _params[3];

    if (features->meltpoint < _params[5]) {
      score -= fabs(features->meltpoint - _params[5]) * _params[4];
    } else if (features->meltpoint > _params[6]) {
      score -= fabs(features->meltpoint - _params[6]) * _params[4];
    }
    return score;
  }
};

} // namespace eternabot

#endif /* defined(__RNAMake__a_basic_test__) */
