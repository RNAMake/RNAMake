//
//  simulate_tectos.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/18/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__simulate_tectos__
#define __RNAMake__simulate_tectos__

#include <stdio.h>

//RNAMake Headers
#include "base/application.hpp"
#include "thermo_fluctuation/thermo_fluc_simulation_devel.h"
#include "motif_data_structure/motif_state_ensemble_tree.h"



class SimulateTectosAppException : public std::runtime_error {
public:
    SimulateTectosAppException(
        String const & message):
    std::runtime_error(message)
    {}
};


String
remove_Us(String const &);


class SimulateTectosApp : public base::Application {
public:
    SimulateTectosApp();
    
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
    
    //old code for comparision
    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_old(
        String const &,
        String const &,
        String const &,
        String const &);
    
    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_new_receptor(
        String const &,
        String const &,
        String const &,
        String const &);
    
    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_old_full_seq(
        String const &,
        String const &,
        String const &,
        String const &);
    
    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_old_reverse(
        String const &,
        String const &,
        String const &,
        String const &);
    
    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_old_coorigin(
        String const &,
        String const &,
        String const &,
        String const &);

    motif_data_structure::MotifStateEnsembleTreeOP
    get_mset_old_gaaa_coorigin(
            String const &,
            String const &,
            String const &,
            String const &);

    thermo_fluctuation::ThermoFlucSimulationLoggerOP
    _get_logger(
            String const &);

    thermo_fluctuation::ThermoFlucScorerOP
    _get_scorer(
            String const &);

    motif::MotifOPs
    get_motifs_from_seq_and_ss(
        String const &,
        String const &);
    
    
private:
    thermo_fluctuation::ThermoFlucSimulationDevel tfs_;
    Strings motif_names_;
    Strings ggaa_ttr_end_names_;
    Strings gaaa_ttr_end_names_;
    
    
    
};

#endif /* defined(__RNAMake__simulate_tectos__) */
