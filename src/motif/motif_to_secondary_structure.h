//
//  motif_to_secondary_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_to_secondary_structure__
#define __RNAMake__motif_to_secondary_structure__

#include <stdio.h>
#include <map>
#include <queue>

#include "secondary_structure/chain.h"
#include "motif/motif.h"

namespace motif {

class MotiftoSecondaryStructure {
public:
    MotiftoSecondaryStructure() :
            chains_(structure::ChainOPs()),
            open_chains_(std::queue<structure::ChainOP>()),
            seen_res_(std::map<util::Uuid, int, util::UuidCompare>()),
            seen_bp_(std::map<util::Uuid, structure::BasepairOP, util::UuidCompare>()) {}

public:

    inline
    void
    reset() {
        chains_ = structure::ChainOPs();
        open_chains_ = std::queue<structure::ChainOP>();
        seen_res_ = std::map<util::Uuid, int, util::UuidCompare>();
        seen_bp_ = std::map<util::Uuid, structure::BasepairOP, util::UuidCompare>();
    }

    secondary_structure::RNAStructureOP
    to_secondary_structure(
            structure::RNAStructureOP const &);

private:

    structure::ChainOP
    _get_next_chain(
            structure::RNAStructureOP const &);

    secondary_structure::RNAStructureOP
    _setup_basepairs_and_ends(
            secondary_structure::StructureOP &,
            structure::RNAStructureOP const &);

private:
    structure::ChainOPs chains_;
    std::queue<structure::ChainOP> open_chains_;
    std::map<util::Uuid, int, util::UuidCompare> seen_res_;
    std::map<util::Uuid, structure::BasepairOP, util::UuidCompare> seen_bp_;


};

}

#endif /* defined(__RNAMake__motif_to_secondary_structure__) */
