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
#include "structure/residue_type_set.h"
#include "structure/residue.h"
#include "structure/basepair.h"
#include "structure/structure.h"
#include "motif/motif_type.h"

class Motif {
public:
    Motif():
    beads_(Beads()),
    score_(0),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
    mdir_(String()),
    name_(String()),
    cached_rotations_(Matrices()),
    structure_ ( StructureOP() )
    {}
    
    Motif(String const &);
    
    Motif(String const &,
          ResidueTypeSet const &);
    
    
    Motif(StructureOP const &,
          BasepairOPs const &);
    
    Motif
    copy();
    
    ~Motif() {}

public:
    
    inline
    void
    _cache_basepair_frames() {
        cached_rotations_ = Matrices(basepairs_.size());
        int i = 0;
        for (auto & bp : basepairs_ ) {
            cached_rotations_[i] = bp->r();
            i++;
        }
    }
    
    String
    sequence();
    
    String
    secondary_structure();
    
    BasepairOPs
    get_basepair(Uuid const &);
    
    BasepairOPs
    get_basepair(ResidueOP,
                 ResidueOP);
    
    BasepairOPs
    get_basepair(Uuid const &,
                 Uuid const &);
    
    inline
    ResidueOP const
    get_residue(int num,
                String const & chain_id,
                String const & i_code) {
        return structure_->get_residue(num, chain_id, i_code);
    }
    
    inline
    ResidueOP const
    get_residue(Uuid const & uuid) {
        return structure_->get_residue(uuid);
    }
    
    Beads const &
    get_beads(BasepairOPs const &);
    
    Beads const &
    get_beads(BasepairOP const &);
    
    Beads const &
    get_beads() {
        ResidueOPs res;
        beads_ = structure_->get_beads(res);
        return beads_;
    }
    
    inline
    AtomOPs const
    atoms() { return structure_->atoms(); }
    
    inline
    ResidueOPs const
    residues() { return structure_->residues(); }
    
    ChainOPs const &
    chains() { return structure_->chains(); }
    
public:
    
    inline
    BasepairOPs const &
    ends() const { return ends_; }
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    BasepairOPs const &
    basepairs() const { return basepairs_; }
    
    inline
    Beads const &
    beads() const { return beads_; }
    
    inline
    MotifType const &
    mtype() { return mtype_; }
    
    inline
    StructureOP const &
    structure() { return structure_; }
    
public: // setters
    
    inline
    void
    structure(StructureOP const & nstructure) { structure_ = nstructure; }
    
public:
    
    void
    setup_basepair_ends();
    
    inline
    void
    transform(Transform const & t) {
        Matrix r_T = t.rotation();
        Matrix transformed;
        r_T.transpose();
        for (auto & bp : basepairs_) {
            dot(bp->r(), r_T, transformed);
            bp->r(transformed);
        }
        structure_->transform(t);
    }
    
    inline
    void
    move(Point const & p) { structure_->move(p); }

    inline
    void
    reset() {
        int i = 0;
        for (auto const & bp : basepairs_) {
            bp->r(cached_rotations_[i]);
            i++;
        }
        
        for (auto const & end : ends_) { end->flip(0); }

        structure_->_restore_coords();
        beads_ = Beads();
    }
    
    String const
    to_str();
    
    String const
    to_pdb_str();
    
    void
    to_pdb(String const);
    
    
    
protected:
    Beads beads_;
    float score_;
    MotifType mtype_;
    BasepairOPs basepairs_, ends_;
    String mdir_, name_;
    Matrices cached_rotations_;
    StructureOP structure_;
};


typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP>   MotifOPs;

void
align_motif(BasepairOP const &,
            BasepairOP const &,
            MotifOP const &);

Motif
ref_motif();


#endif /* defined(__RNAMake__motif__) */
