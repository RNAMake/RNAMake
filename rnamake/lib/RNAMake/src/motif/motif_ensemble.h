//
//  motif_ensemble.h
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_ensemble__
#define __RNAMake__motif_ensemble__

#include <stdio.h>

//RNAMake Headers
#include "motif/motif.h"

struct MotifEnsembleMember {
    MotifEnsembleMember(
        MotifOP const & nmotif,
        float const & nenergy):
    motif(nmotif),
    energy(nenergy)
    {}
    
    MotifOP motif;
    float energy;
};

typedef std::shared_ptr<MotifEnsembleMember> MotifEnsembleMemberOP;
typedef std::vector<MotifEnsembleMemberOP>   MotifEnsembleMemberOPs;

class MotifEnsemble {
public:
    MotifEnsemble():
    id_(""),
    block_end_add_(0),
    members_(MotifEnsembleMemberOPs())
    {}
    
private:
    String id_;
    int block_end_add_;
    MotifEnsembleMemberOPs members_;
    
    
};



#endif /* defined(__RNAMake__motif_ensemble__) */
