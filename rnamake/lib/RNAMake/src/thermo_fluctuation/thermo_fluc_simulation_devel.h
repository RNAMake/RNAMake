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
#include "base/option.h"
#include "thermo_fluctuation/thermo_fluc_sampler.h"
#include "thermo_fluctuation/thermo_fluc_scorer.h"


class ThermoFlucSimulationException : public std::runtime_error {
public:
    ThermoFlucSimulationException(
        String const & message) :
    std::runtime_error(message)
    {}
};


class ThermoFlucSimulationDevel : public OptionClass {
public:
    ThermoFlucSimulationDevel() {
        scorer_ = std::make_shared<FrameScorerDevel>(FrameScorerDevel());
        sampler_ = ThermoFlucSampler();
        check_nodes_1_ = {  };
        check_nodes_2_ = {  };
        record_state_ = 0;
        record_all_ = 0;
        steric_radius_ = 2.2;
        setup_options();
        
    }
    
    ~ThermoFlucSimulationDevel() {}

private:
    inline
    int
    _check_sterics() {
        clash_ = 0;
        for(auto const & i : check_nodes_1_) {
            for(auto const & j : check_nodes_2_) {
                for(auto const & b2 : sampler_.mst()->get_node(i)->data()->cur_state->beads()) {
                    for(auto const & b1 : sampler_.mst()->get_node(j)->data()->cur_state->beads()) {
                        if(b1.distance(b2) < steric_radius_) { clash_ = 1; }
                    }
                    
                    if(clash_) { break; }
                }
                if(clash_) { break; }
            }
            if(clash_) { break; }
        }
        return clash_;
    }
    
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
    
   
public: //option wrappers
    
    inline
    Options &
    options() { return options_; }
    
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
    
    void
    update_var_options();
    
    
protected:
    void
    setup_options();
    
    
    
private:
    ThermoFlucScorerOP scorer_;
    ThermoFlucSampler sampler_;
    BasepairStateOP end_state_1_, end_state_2_;
    Options options_;
    int ni1_, ni2_, ei1_, ei2_;
    int clash_;
    float score_;
    Ints check_nodes_1_, check_nodes_2_;
    //option vars
    String record_file_, record_all_file_;
    float temperature_, cutoff_, steric_radius_;
    int steps_, record_state_, record_all_;
    bool record_;
    int bound_pdb_;
};

#endif /* defined(__RNAMake__thermo_fluc_simulation_devel__) */
