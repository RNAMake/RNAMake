//
//  motif.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif__
#define __RNAMake__motif__

#include <stdio.h>


//RNAMake Headers
#include "base/types.h"
#include "math/xyz_matrix.h"
#include "math/transform.h"
#include "util/uuid.h"
#include "secondary_structure/motif.h"
#include "structure/residue_type_set.h"
#include "structure/residue.h"
#include "structure/basepair.h"
#include "structure/structure.h"
#include "structure/rna_structure.h"
#include "util/motif_type.h"
#include "motif/motif_state.h"

namespace motif {

class Motif : public structure::RNAStructure {
public:
    Motif() :
            structure::RNAStructure(),
            id_(util::Uuid()),
            block_end_add_(0) {}

    Motif(
            structure::StructureOP const & structure,
            structure::BasepairOPs const & basepairs,
            structure::BasepairOPs const & ends) :
            structure::RNAStructure(structure, basepairs, ends),
            id_(util::Uuid()),
            block_end_add_(0) {}

    Motif(
            String const &,
            structure::ResidueTypeSet const &);


    Motif(
            Motif const & m);

    Motif(
            structure::RNAStructure const & rs) :
            structure::RNAStructure(rs) {
        id_ = util::Uuid();
        block_end_add_ = -1;
        secondary_structure_ = std::make_shared<secondary_structure::Motif>();
        mtype_ = util::MotifType::UNKNOWN;
    }

    ~Motif() {}

public:


    inline
    void
    transform(
            math::Transform const & t) {

        math::Matrix r_T = t.rotation();
        math::Matrix transformed;
        r_T.transpose();
        for (auto & bp : basepairs_) {
            dot(bp->r(), r_T, transformed);
            bp->r(transformed);
        }
        structure_->transform(t);
    }

    inline
    void
    move(
            math::Point const & p) { structure_->move(p); }

    String const
    to_str();


    MotifStateOP
    get_state();


    void
    new_res_uuids();

    void
    copy_uuids_from_motif(
            Motif const &);

public: //wrappers from secondary structure

    String
    sequence() { return secondary_structure_->sequence(); }

    String
    dot_bracket() { return secondary_structure_->dot_bracket(); }

public: //getters

    inline
    util::MotifType const &
    mtype() { return mtype_; }

    inline
    secondary_structure::MotifOP const &
    secondary_structure() { return secondary_structure_; }

    inline
    int const &
    block_end_add() { return block_end_add_; }

    inline
    util::Uuid const &
    id() { return id_; }

    inline
    String
    end_name(int i) { return ends_[i]->name(); }


public: // setters

    inline
    void
    id(util::Uuid const & nid) { id_ = nid; }

    inline
    void
    mtype(util::MotifType const & mtype) { mtype_ = mtype; }

    inline
    void
    secondary_structure(
            secondary_structure::MotifOP const & ss) {
        secondary_structure_ = ss;
    }

    inline
    void
    structure(
            structure::StructureOP const & s) {
        structure_ = s;
    }

    inline
    void
    block_end_add(
            int nblock_end_add) {
        block_end_add_ = nblock_end_add;
    }


private:
    util::MotifType mtype_;
    secondary_structure::MotifOP secondary_structure_;
    int block_end_add_;
    util::Uuid id_;
};


typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP>   MotifOPs;

void
align_motif(
        structure::BasepairStateOP const &,
        structure::BasepairOP const &,
        MotifOP &);

MotifOP
get_aligned_motif(
        structure::BasepairOP const &,
        structure::BasepairOP const &,
        MotifOP const &);


Motif
ref_motif();

MotifOP
file_to_motif(
        String const &);

inline
int
clash_between_motifs(
        MotifOP const & m1,
        MotifOP const & m2,
        double clash_radius = 2.7) {

    for (auto const & b1 : m1->beads()) {
        if (b1.btype() == structure::BeadType::PHOS) { continue; }
        for (auto const & b2 : m2->beads()) {
            if (b2.btype() == structure::BeadType::PHOS) { continue; }
            if (b1.distance(b2) < clash_radius) { return 1; }
        }
    }
    return 0;
}

}


#endif /* defined(__RNAMake__motif__) */
