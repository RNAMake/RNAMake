//
//  pose.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 2/17/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure/chain.h"
#include "motif/pose.h"

Pose::Pose(MotifOP const & m) {
    structure_ = m->structure();
    basepairs_ = m->basepairs();
    ends_ = m->ends();
    designable_ = std::map<Uuid, int, UuidCompare> ();
}

Pose::Pose(
    StructureOP const & structure,
    BasepairOPs const & basepairs) {
    structure_ = structure;
    basepairs_ = basepairs;
    designable_ = std::map<Uuid, int, UuidCompare> ();
}

String
Pose::designable_sequence() {
    String seq, s;
    for (auto const & c : chains()) {
        for (auto const & r : c->residues()) {
            BasepairOPs bps = get_basepair(r->uuid());
            s = r->short_name();
            for (auto const & bp : bps) {
                if(designable_.find(bp->uuid()) != designable_.end() ){
                    s = "N"; break;
                }
            }
            seq += s;
        }
        seq += "&";
    }
    
    return seq.substr(0, seq.length()-1);
}



void
Pose::set_bp_designable(BasepairOP const & bp) {
    designable_[bp->uuid()] = 1;
}