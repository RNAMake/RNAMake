//
//  monte_carlo.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__monte_carlo__
#define __RNAMake__monte_carlo__

#include <math.h>
#include <stdio.h>

// RNAMake Headers
#include "util/random_number_generator.h"

namespace util {

class MonteCarlo {
public:
  inline MonteCarlo(float temperature = 1.0f)
      : _temperature(temperature), _rng(RandomNumberGenerator()) {}

  ~MonteCarlo() {}

  inline int accept(float current, float next) {

    if (next < current) {
      return 1;
    }

    _score = exp((current - next) / _temperature);
    if (_rng.rand() < _score) {
      return 1;
    }

    return 0;
  }

public:
  inline void set_temperature(float new_temp) { _temperature = new_temp; }

  inline float get_temperature() { return _temperature; }

  inline void scale_temperature(float scale) { _temperature *= scale; }

private:
  float _temperature;
  float _score;
  RandomNumberGenerator _rng;
};

} // namespace util

#endif /* defined(__RNAMake__monte_carlo__) */
