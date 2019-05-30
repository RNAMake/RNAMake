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
#include "motif_data_structure/motif_state_ensemble_tree.h"


namespace thermo_fluctuation {

class ThermoFlucSamplerException : public std::runtime_error {
public:
    ThermoFlucSamplerException(
            String const & message) :
            std::runtime_error(message) {}
};

class ThermoFlucSampler {
public:
    ThermoFlucSampler():
            temperature_(298.15 * 1.3806488e-1),
            rng_(util::RandomNumberGenerator()),
            mc_(util::MonteCarlo(temperature_)) {}

    ~ThermoFlucSampler() {}

public:

    void
    sample(
            int steps = 1000);

    void
    setup(
            motif_data_structure::MotifStateEnsembleTreeOP const &);

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
            motif::MotifStateEnsembleMemberOP const &);


public: // getters

    inline
    float
    temperature() { return temperature_; }

    inline
    motif_data_structure::MotifStateTreeOP
    mst() { return mst_; }

public: // setters
    inline
    void
    temperature(float const & temp) {
        if (temp < 0) {
            throw ThermoFlucSamplerException(
                    "cannot set temperature lower then 0");
        }

        temperature_ = temp * 1.3806488e-1;
        mc_ = util::MonteCarlo(temperature_);
    }

    inline
    void
    randomized_start(bool random_start) {
        randomized_start_ = random_start;
    }

private:
    bool randomized_start_;
    float temperature_;
    util::MonteCarlo mc_;
    util::RandomNumberGenerator rng_;
    motif_data_structure::MotifStateEnsembleTreeOP mset_;
    motif_data_structure::MotifStateTreeOP mst_;
    Ints states_;
    int node_num_, pos_, mem_pos_, accept_;
    float energy_;
    motif_data_structure::MotifStateEnsembleTreeNodeOP mset_node_;
    motif_data_structure::MotifStateTreeNodeOP mst_node_;
    motif::MotifStateEnsembleMemberOP new_mem_;
    //hold last move to undo
    int last_state_pos_, last_num_;
    motif::MotifStateOP last_state_;
};

}

#endif /* defined(__RNAMake__thermo_fluc_sampler__) */
