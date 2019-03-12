//
//  thermo_fluc_sampler.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/6/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__thermo_fluc_sampler__
#define __RNAMake__thermo_fluc_sampler__

#include <stdio.h>

//RNAMake Headers
#include "util/monte_carlo.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"


class ThermoFlucSamplerException : public std::runtime_error {
public:
    ThermoFlucSamplerException(
        String const & message) :
    std::runtime_error(message)
    {}
};

class ThermoFlucSampler {
public:
    ThermoFlucSampler():
    temperature_(298.15),
    rng_(util::RandomNumberGenerator())
    {}
    
    ~ThermoFlucSampler() {}
    
public:
    
    void
    setup(
        MotifStateEnsembleTreeOP const &);
    
    void
    sample(
        int steps=1000);
    
    int
    next();
    
    void
    undo();
    
    void
    to_pdb(
        String fname = "test.pdb",
        int renumber = -1);
    
private:
    
    void
    update(
        int,
        MotifStateEnsembleMemberOP const &);
    

    
public: // getters
    
    inline
    float
    temperature() { return temperature_; }
    
    inline
    MotifStateTreeOP
    mst() { return mst_; }
    
public: // setters
    inline
    void
    temperature(float const & temp) {
        if(temp < 0) {
            throw ThermoFlucSamplerException(
                "cannot set temperature lower then 0");
        }
        
        temperature_ = temp; }
    
private:
    float temperature_;
    util::MonteCarlo mc_;
    util::RandomNumberGenerator rng_;
    MotifStateEnsembleTreeOP mset_;
    MotifStateTreeOP mst_;
    Ints states_;
    int node_num_, pos_, mem_pos_, accept_;
    float energy_;
    MotifStateEnsembleTreeNodeOP mset_node_;
    MotifStateTreeNodeOP mst_node_;
    MotifStateEnsembleMemberOP new_mem_;
    //hold last move to undo
    int last_state_pos_, last_num_;
    MotifStateOP last_state_;
    
    
};

#endif /* defined(__RNAMake__thermo_fluc_sampler__) */
