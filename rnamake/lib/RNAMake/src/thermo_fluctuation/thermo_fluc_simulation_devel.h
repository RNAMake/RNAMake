//
//  thermo_fluc_simulation_devel.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_simulation_devel__
#define __RNAMake__thermo_fluc_simulation_devel__

#include <stdio.h>

#include "base/types.h"
#include "base/base.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"

class ThermoFlucSimulationDevel : public Base {
public:
    ThermoFlucSimulationDevel() {
        scorer_ = std::make_shared<FrameScorerDevel>(FrameScorerDevel());
        sampler_ = ThermoFlucSampler();
        setup_options();
        
    }
    
    ~ThermoFlucSimulationDevel() {}
    
public:
    
    void
    setup(
          MotifStateEnsembleTreeOP const &,
          int,
          int,
          int,
          int);
    
    int
    run();
    
    String
    static_run();
    
protected:
    void
    setup_options();
    
private:
    
    void
    update_var_options();
    
    
private:
    ThermoFlucScorerOP scorer_;
    ThermoFlucSampler sampler_;
    BasepairStateOP end_state_1_, end_state_2_;
    int ni1_, ni2_, ei1_, ei2_;
    float score_;
    //option vars
    String record_file_;
    float temperature_, cutoff_, steric_radius_;
    int steps_, record_;
};


#endif /* defined(__RNAMake__thermo_fluc_simulation_devel__) */
