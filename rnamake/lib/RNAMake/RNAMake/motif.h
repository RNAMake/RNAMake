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
#include "types.h"
#include "motif_type.h"
#include "residue.h"
#include "basepair.h"
#include "structure.h"
#include "residue_type_set.h"
#include "transform.h"
#include "xyzMatrix.h"

class Motif {
public:
    Motif():
    beads_(Beads()),
    score_(0),
    basepairs_(Basepairs()),
    ends_(Basepairs()),
    mdir_(String()),
    name_(String()),
    cached_rotations_(Matrices())
    {}
    
    Motif(String const &,
          ResidueTypeSet const &);
    
    Motif
    copy();
    
    ~Motif() {}

public:
    inline
    void
    _cache_basepair_frames() {
        cached_rotations_ = Matrices(basepairs_.size());
        int i = 0;
        for (auto const & bp : basepairs_ ) {
            cached_rotations_[i] = bp.r();
            i++;
        }
    }
    
    Basepairs
    get_basepair(Uuid const &);
    
    Basepairs
    get_basepair(Residue const &,
                 Residue const &);
    
    Basepairs
    get_basepair(Uuid const &,
                 Uuid const &);
    
    inline
    Residue const
    get_residue(int num,
                String const & chain_id,
                String const & i_code) {
        return structure_.get_residue(num, chain_id, i_code);
    }
    
    inline
    Residue const
    get_residue(Uuid const & uuid) {
        return structure_.get_residue(uuid);
    }
    
    inline
    Residues const
    residues() { return structure_.residues(); }
    
    inline
    Chains const
    chains() { return structure_.chains(); }
    
    inline
    Basepairs const &
    ends() const { return ends_; }
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    Basepairs const &
    basepairs() const { return basepairs_; }
    
    
public:
    
    inline
    void
    transform(Transform const & t) {
        Matrix r_T = t.rotation();
        Matrix transformed;
        r_T.transpose();
        for (auto & bp : basepairs_) {
            dot(bp.r(), r_T, transformed);
            bp.r(transformed);
        }
        structure_.transform(t);
    }
    
    inline
    void
    move(Point const & p) { structure_.move(p); }

    String const
    to_str();
    
    String const
    to_pdb_str();
    
    void
    to_pdb(String const);
    
    
    
private:
    Beads beads_;
    float score_;
    MotifType mtype_;
    Basepairs basepairs_, ends_;
    String mdir_, name_;
    Matrices cached_rotations_;
    Structure structure_;
};

void
align_motif(Basepair const &,
            Basepair const &,
            Motif &);

#endif /* defined(__RNAMake__motif__) */
