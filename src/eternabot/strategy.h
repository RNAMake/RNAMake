//
//  strategy.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__strategy__
#define __RNAMake__strategy__

#include <stdio.h>

#include "eternabot/feature_generator.h"

namespace eternabot {

class Strategy {
public:
  Strategy() {}

  virtual ~Strategy() {}

public:
  virtual float score(FeaturesOP const &) = 0;

public:
  inline float mean() const { return _mean; }

  inline float stdev() const { return _stdev; }

protected:
  float _mean, _stdev;
  Reals _params;
};

typedef std::shared_ptr<Strategy> StrategyOP;
typedef std::vector<StrategyOP> StrategyOPs;

} // namespace eternabot

#endif /* defined(__RNAMake__strategy__) */
