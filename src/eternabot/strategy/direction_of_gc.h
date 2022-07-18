//
//  direction_of_gc.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/4/16.
//  Copyright (c) 2016 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__direction_of_gc__
#define __RNAMake__direction_of_gc__

#include <stdio.h>

#include <eternabot/strategy.h>

namespace eternabot {

class DirectionofGCPairsinMultiLoops : public Strategy {
public:
  DirectionofGCPairsinMultiLoops() {
    _mean = 85.2869664088;
    _stdev = 26.9535204308;
  }

  ~DirectionofGCPairsinMultiLoops() {}

  float score(FeaturesOP const &features) { return 100; }
};

} // namespace eternabot

#endif /* defined(__RNAMake__direction_of_gc__) */
