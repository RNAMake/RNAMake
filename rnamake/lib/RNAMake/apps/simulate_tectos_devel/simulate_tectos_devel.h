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
#include "base/cl_option.h"
#include "motif_data_structures/motif_state_ensemble_tree.h"


struct MotifInfo {
    String name, end_id, end_name;
};

typedef std::vector<MotifInfo> MotifInfos;



class SimulateTectos {
public:
    SimulateTectos(
        Options &);
    
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
    
    MotifStateEnsembleTreeOP
    get_mset_old_full_seq(
        String const &,
        String const &,
        String const &,
        String const &);
    
    MotifStateEnsembleTreeOP
    get_mset_old_reverse(
        String const &,
        String const &,
        String const &,
        String const &);
    
    MotifStateEnsembleTreeOP
    get_mset_old_coorigin(
        String const &,
        String const &,
        String const &,
        String const &);
    
    MotifInfos
    get_motifs_from_seq_and_ss(
        String const &,
        String const &);
    
    
    
    
};

#endif /* defined(__RNAMake__simulate_tectos__) */
