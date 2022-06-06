//
//  motif_scorer.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/16/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_scorer__
#define __RNAMake__motif_scorer__

#include <map>
#include <stdio.h>

//RNAMake Headers
#include "base/types.hpp"
#include "structure/basepair.h"
#include "motif/motif.h"

namespace motif {

class MotifScorer {
public:
    MotifScorer();

    ~MotifScorer() {}

public:
    float
    score(
            MotifOP const &);

private:

    float
    _score_cWW_bp(
            structure::BasepairOP const &);

    void
    _bp_reference_energy_table();

private:
    StringFloatMap bp_ref_energy_;
    float unpaired_pentalty_;
};

}

#endif /* defined(__RNAMake__motif_scorer__) */
