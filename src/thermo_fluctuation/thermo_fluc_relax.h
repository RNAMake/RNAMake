//
// Created by Joseph Yesselman on 4/25/18.
//

#ifndef TEST_THERMO_FLUC_RELAX_H
#define TEST_THERMO_FLUC_RELAX_H

#include "thermo_fluctuation/thermo_fluc_sampler.h"

namespace thermo_fluctuation {

class ThermoFlucRelax {
public:
  ThermoFlucRelax() : sampler_(ThermoFlucSampler()) {
    sampler_.temperature(1000);
  }

  ~ThermoFlucRelax() {}

public:
  void run_with_graph();

private:
  ThermoFlucSampler sampler_;
};

} // namespace thermo_fluctuation

#endif // TEST_THERMO_FLUC_RELAX_H
