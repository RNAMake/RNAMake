//
// Created by Joseph Yesselman on 8/30/17.
//

#ifndef TEST_GENERAL_HELIX_SAMPLER_H
#define TEST_GENERAL_HELIX_SAMPLER_H

#include <stdio.h>

#include "base/application.hpp"
#include "thermo_fluctuation/thermo_fluc_simulation.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"

class GeneralHelixSampler : public Application {
public:
    GeneralHelixSampler();

    ~GeneralHelixSampler() {}

public:

    void
    setup_options();

    void
    parse_command_line(
            int,
            const char **);

    void
    run();

private:
    MotifOPs
    get_motifs_from_seq_and_ss(
            String const &,
            String const &);

private:

    ThermoFlucSimulation tfs_;

};





#endif //TEST_GENERAL_HELIX_SAMPLER_H
