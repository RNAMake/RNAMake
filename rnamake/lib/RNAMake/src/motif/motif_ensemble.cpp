//
//  motif_ensemble.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_ensemble.h"

namespace motif {

MotifStateEnsembleOP
MotifEnsemble::get_state() {
    auto motif_states = MotifStateOPs();
    auto energies = Floats();

    for (auto const & mem : members_) {
        motif_states.push_back(mem->motif->get_state());
        energies.push_back(mem->energy);
    }

    return std::make_shared<MotifStateEnsemble>(motif_states, energies);

}

}