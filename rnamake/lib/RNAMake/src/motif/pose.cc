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
    /*structure_ = m->structure();
    basepairs_ = m->basepairs();
    secondary_structure_ = m->secondary_structure();
    ends_ = m->ends();
    designable_ = std::map<util::Uuid, int, util::UuidCompare> ();
    path_ = m->path();*/
}

Pose::Pose(
    structure::StructureOP const & structure,
    structure::BasepairOPs const & basepairs) {
    structure_ = structure;
    basepairs_ = basepairs;
    designable_ = std::map<util::Uuid, int, util::UuidCompare> ();
}

MotifOPs const &
Pose::motifs(
    util::MotifType const & mtype) {
    
    if(motifs_.find(mtype) == motifs_.end()) {
        motifs_[mtype] = MotifOPs();
    }
    return motifs_[mtype];

}

String
Pose::designable_sequence() {
    String seq, s;
    for (auto const & c : chains()) {
        for (auto const & r : c->residues()) {
            structure::BasepairOPs bps = get_basepair(r->uuid());
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
Pose::set_bp_designable(structure::BasepairOP const & bp) {
    designable_[bp->uuid()] = 1;
}