//
//  pose.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/7/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake_sec_pose__
#define __RNAMake_sec_pose__

#include <stdio.h>

// RNAMake Headers
#include "secondary_structure/rna_structure.h"
#include "secondary_structure/motif.h"

namespace secondary_structure {

class Pose : public RNAStructure {
public:
    Pose() :
            RNAStructure() {}

    Pose(
            StructureOP const & structure,
            BasepairOPs const & basepairs,
            BasepairOPs const & ends) :
            RNAStructure(structure, basepairs, ends) {}

    Pose(
            StructureOP const & structure,
            BasepairOPs const & basepairs,
            BasepairOPs const & ends,
            MotifOPs const & motifs) :
            RNAStructure(structure, basepairs, ends),
            motifs_(motifs) {}

    Pose(
            RNAStructureOP const & rs,
            MotifOPs const & motifs) :
            motifs_(motifs) {
        this->structure_ = rs->structure();
        this->basepairs_ = rs->basepairs();
        this->ends_ = rs->ends();
        this->end_ids_ = rs->end_ids();
    }

public:
    MotifOPs const &
    helices() {
        if (helices_.size() == 0) {
            _build_helices();
        }
        return helices_;
    }

    MotifOPs const &
    motifs() { return motifs_; }

    MotifOP
    motif(
            util::Uuid const & uuid) {
        for (auto const & m : motifs_) {
            if (m->id() == uuid) { return m; }
        }
        return nullptr;
    }

    void
    replace_sequence(
            String const &);

    void
    update_motif(
            util::Uuid const &);


private:
    void
    _build_helices();

private:
    MotifOPs motifs_;
    MotifOPs helices_;


};

typedef std::shared_ptr<Pose> PoseOP;

}
#endif /* defined(__RNAMake__pose__) */
