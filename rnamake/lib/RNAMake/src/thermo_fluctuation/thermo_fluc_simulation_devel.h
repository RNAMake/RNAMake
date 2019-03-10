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

class ThermoFlucSimulationLogger {
public:
    ThermoFlucSimulationLogger(
            String const & fname = "test.out"):
            outputed_header_(false){
        out_.open(fname);
    }

    ~ThermoFlucSimulationLogger() {
        out_.close();
    }

public:
    virtual
    void
    log(
            MotifStateTreeOP const & mst,
            float score) {

        auto end_state_1 = mst->get_node(ni1_)->data()->get_end_state(ei1_);
        auto end_state_2 = mst->get_node(ni2_)->data()->get_end_state(ei2_);

        out_ << vector_to_str(end_state_1->d()) << "," << matrix_to_str(end_state_1->r()) << ",";
        out_ << vector_to_str(end_state_2->d()) << "," << matrix_to_str(end_state_2->r()) << ","  << score;
        out_ << std::endl;
    }

    virtual
    void
    setup(
            MotifStateTreeOP const & mst,
            int ni1,
            int ni2,
            int ei1,
            int ei2) {
        ni1_ = ni1;
        ni2_ = ni2;
        ei1_ = ei1;
        ei2_ = ei2;
        if(! outputed_header_) {
            _output_header(mst);
        }
    }

    virtual
    void
    finalize() {

    }

protected:
    virtual
    void
    _output_header(
            MotifStateTreeOP const & mst) {
        out_ << "d_target,r_target,d_variable,r_variable,score" << std::endl;
        outputed_header_ = true;
    }

protected:
    std::ofstream out_;
    int ni1_, ni2_, ei1_, ei2_;
    bool outputed_header_;


};

typedef std::shared_ptr<ThermoFlucSimulationLogger> ThermoFlucSimulationLoggerOP;

class ThermoFlucSimulationDevel : public base::OptionClass {
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

public:
    inline
    void
    set_logger(
            ThermoFlucSimulationLoggerOP logger) {
        logger_ = logger;
    }

    inline
    void
    set_scorer(
            ThermoFlucScorerOP scorer) {
        scorer_ = scorer;
    }
   
public: //option wrappers
    
    inline
    base::Options &
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
    ThermoFlucSimulationLoggerOP logger_;
    ThermoFlucSampler sampler_;
    BasepairStateOP end_state_1_, end_state_2_;
    base::Options options_;
    int ni1_, ni2_, ei1_, ei2_;
    int clash_;
    float score_;
    Ints check_nodes_1_, check_nodes_2_;
    //option vars
    String record_file_, record_all_file_;
    float temperature_, cutoff_, steric_radius_;
    int steps_, record_state_, record_all_;
    bool record_;
    bool record_only_bound_, record_only_unbound_;
    bool dump_state_;
    bool dump_pdbs_;
};

#endif /* defined(__RNAMake__thermo_fluc_simulation_devel__) */
