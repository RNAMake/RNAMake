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
#include "base/base.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"

String
bp_name_from_sequence(
    String const &);


class SimulateTectos {
public:
    SimulateTectos(
        String const &,
        String const &,
        String const &,
        String const &);
    
public:
    
    void
    run();

private:
    
    //old code for comparision
    MotifStateEnsembleTreeOP
    get_mset_old(
        String const &,
        String const &,
        String const &,
        String const &);
    
    Strings
    get_motifs_from_seq_and_ss(
        String const &,
        String const &);
    
    
    
    
};

#endif /* defined(__RNAMake__simulate_tectos__) */
