//
//  motif_state_ensemble.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 9/3/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "motif/motif_state_ensemble.h"

namespace motif {

MotifStateEnsemble::MotifStateEnsemble(
        String const & s) :
        id_(""),
        block_end_add_(0),
        members_(MotifStateEnsembleMemberOPs()) {

    auto spl = base::split_str_by_delimiter(s, "{");
    id_ = spl[0];
    block_end_add_ = std::stoi(spl[1]);
    for (int i = 2; i < spl.size(); i++) {
        auto spl2 = base::split_str_by_delimiter(spl[i], "#");
        auto ms = std::make_shared<MotifState>(spl2[0]);
        auto energy = std::stof(spl2[1]);
        auto mem = std::make_shared<MotifStateEnsembleMember>(ms, energy);
        members_.push_back(mem);
    }
}

MotifStateEnsemble::MotifStateEnsemble(
        MotifStateOPs const & motif_states,
        Floats const & energies) {

    if (motif_states.size() != energies.size()) {
        throw MotifStateEnsembleException("must supply the same number of motif_states and energies");
    }

    if (motif_states.size() == 0) {
        throw MotifStateEnsembleException("must supply at least one motif state to ensemble");
    }

    for (int i = 0; i < motif_states.size(); i++) {
        auto member = std::make_shared<MotifStateEnsembleMember>(motif_states[i], energies[i]);
        members_.push_back(member);
    }

    std::sort(members_.begin(), members_.end(), MotifStateEnsembleMember_LessThanKey());
    block_end_add_ = members_[0]->motif_state->block_end_add();
    id_ = members_[0]->motif_state->end_ids()[0];

}


String
MotifStateEnsemble::to_str() {
    String s = id_ + "{" + std::to_string(block_end_add_) + "{";
    for (auto const & mem : members_) {
        s += mem->to_str() + "{";
    }
    return s;
}

}

// non-member functions     ////////////////////////////////////////////////////////////////////////
