//
//  thermo_fluc_simulation.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/8/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_simulation__
#define __RNAMake__thermo_fluc_simulation__

#include <stdio.h>

#include "base/types.h"
#include "base/base.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"


class ThermoFlucSimulation : public Base {
public:
    ThermoFlucSimulation() {
        scorer_ = std::make_shared<ThermoFlucScorer>(FrameScorer());
        sampler_ = ThermoFlucSampler();
        setup_options();
        
    }

    ~ThermoFlucSimulation() {}
    
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
    
protected:
    void
    setup_options();
    
private:
    
    void
    update_var_options();
    
    
private:
    ThermoFlucScorerOP scorer_;
    ThermoFlucSampler sampler_;
    MotifStateTreeNodeOP n1_, n2_;
    int ei1_, ei2_;
    //option vars
    float temperature_, cutoff_;
    int steps_, record_;
};


#endif /* defined(__RNAMake__thermo_fluc_simulation__) */
