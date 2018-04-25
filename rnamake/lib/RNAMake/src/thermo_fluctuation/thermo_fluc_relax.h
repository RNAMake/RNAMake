//
// Created by Joseph Yesselman on 4/25/18.
//

#ifndef TEST_THERMO_FLUC_RELAX_H
#define TEST_THERMO_FLUC_RELAX_H

#include "thermo_fluctuation/thermo_fluc_sampler.h"

class ThermoFlucRelax {
public:
    ThermoFlucRelax():
            sampler_(ThermoFlucSampler()){
        sampler_.temperature(1000);

    }

    ~ThermoFlucRelax() {}

public:
    void
    run_with_graph();

private:
    ThermoFlucSampler sampler_;
};


#endif //TEST_THERMO_FLUC_RELAX_H
