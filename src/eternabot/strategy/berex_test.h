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

#include <eternabot/strategy.h>

namespace eternabot {

class BerexTest : public Strategy {
public:
  BerexTest() {
    _params = Reals(12);
    _params[0] = 0.233331320365;
    _params[1] = 67.1928215814;
    _params[2] = 0.0846125555161;
    _params[3] = 97.8926617326;
    _params[4] = 0.124731583984;
    _params[5] = 76.3152913746;
    _params[6] = -64.8000139266;
    _params[7] = -25.1005947191;
    _params[8] = 1.2097925221;
    _params[9] = 33.3729652592;
    _params[10] = 171.705337159;
    _params[11] = 1.28554296532;
    _mean = 84.0125821249;
    _stdev = 8.91633847502;
  }

  ~BerexTest() {}

public:
  float score(FeaturesOP const &features) {

    float score = 100;
    score -=
        fabsf(features->g_count / features->length - _params[0]) * _params[1];
    score -=
        fabsf(features->u_count / features->length - _params[2]) * _params[3];
    score -=
        fabsf(features->c_count / features->length - _params[4]) * _params[5];

    if (features->fe < _params[6]) {
      score -= fabsf(features->fe - _params[6]) * _params[8];
    } else if (features->fe > _params[7]) {
      score -= fabsf(features->fe - _params[7]) * _params[8];
    }

    if (features->meltpoint < _params[9]) {
      score -= fabsf(features->meltpoint - _params[9]) * _params[11];
    } else if (features->meltpoint > _params[10]) {
      score -= fabsf(features->meltpoint - _params[10]) * _params[11];
    }

    return score;
  }
};

} // namespace eternabot

#endif /* defined(__RNAMake__berex_test__) */
