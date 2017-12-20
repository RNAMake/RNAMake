//
// Created by Joseph Yesselman on 12/18/17.
//

#ifndef TEST_GENERAL_HELIX_SAMPLER_H
#define TEST_GENERAL_HELIX_SAMPLER_H


#include <stdio.h>

//RNAMake Headers
#include "base/application.hpp"
#include "thermo_fluctuation/thermo_fluc_simulation.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"

class GeneralHelixSamplerException : public std::runtime_error {
public:
    GeneralHelixSamplerException(
            String const & message):
            std::runtime_error(message) {}
};


String
remove_Us(String const &);


class GeneralHelixSampler : public Application {
public:
    GeneralHelixSampler();

public: // application setups functions
    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

public:

    void
    run();

private:
    String
    _generate_structure(
            String const &);

    MotifOPs
    get_motifs_from_seq_and_ss(
            String const &,
            String const &);

private:
    ThermoFlucSimulation tfs_;



};



#endif //TEST_GENERAL_HELIX_SAMPLER_H
