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
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"


class ThermoFlucSimulation {
public:
    ThermoFlucSimulation() {
        scorer_ = std::make_shared<FrameScorer>(FrameScorer());
        sampler_ = ThermoFlucSampler();
        options_ = Options("ThermoFlucSimulationOptions");
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
    
public: //option wrappers
    inline
    float
    get_int_option(String const & name) { return options_.get_int(name); }
    
    inline
    float
    get_float_option(String const & name) { return options_.get_float(name); }
    
    inline
    String
    get_string_option(String const & name) { return options_.get_string(name); }
    
    inline
    bool
    get_bool_option(String const & name) { return options_.get_bool(name); }
    
    
    template<typename T>
    void
    set_option_value(
        String const & name,
        T const & val) {
        options_.set_value(name, val);
        update_var_options();
    }
    
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
    Options options_;
    int ni1_, ni2_, ei1_, ei2_;
    float score_;
    //option vars
    bool record_;
    float temperature_, cutoff_;
    int steps_;
};


#endif /* defined(__RNAMake__thermo_fluc_simulation__) */
