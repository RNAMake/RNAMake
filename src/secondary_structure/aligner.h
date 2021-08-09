//
// Created by Joseph Yesselman on 11/4/17.
//

#ifndef RNAMAKE_NEW_SS_ALIGNER_H
#define RNAMAKE_NEW_SS_ALIGNER_H

#include <primitives/aligner.h>
#include <secondary_structure/segment.h>

namespace secondary_structure {

class Aligner : public primitives::Aligner<Segment, Basepair> {
public:
    Aligner() {}

    ~Aligner() {}

public:
    void
    align(
            Basepair const & ref_bp,
            Segment & seg) const { }

    SegmentOP
    get_aligned(
            Basepair const & ref_bp,
            Segment const & seg) const {
        return std::make_shared<Segment>(seg);
    }
};


}

#endif //RNAMAKE_NEW_ALIGNER_H
