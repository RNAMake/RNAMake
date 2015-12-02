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
#include "secondary_structure/secondary_structure.h"
#include "structure/residue_type_set.h"
#include "structure/residue.h"
#include "structure/basepair.h"
#include "structure/structure.h"
#include "util/motif_type.h"
#include "motif/motif_state.h"

class Motif {
public:
    Motif():
    beads_(Beads()),
    score_(0),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
    path_(String()),
    name_(String()),
    end_ids_(Strings()),
    structure_ ( StructureOP() )
    {}
    
    Motif(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    beads_(Beads()),
    score_(0),
    basepairs_(basepairs),
    ends_(ends),
    path_(String()),
    name_(String()),
    end_ids_(Strings()),
    structure_ (structure)
    {}
        
    Motif(String const &,
          ResidueTypeSet const &);
    
    
    Motif
    copy();
    
    ~Motif() {}

public:

    int
    end_index(
        BasepairOP const &);
    
    BasepairOPs 
    get_basepair(
        String const &);
    
    BasepairOPs
    get_basepair(
        Uuid const &);
    
    BasepairOPs
    get_basepair(
        ResidueOP const &,
        ResidueOP const &);
    
    BasepairOPs
    get_basepair(
        Uuid const &,
        Uuid const &);
    
    Beads const &
    get_beads(
        BasepairOPs const &);
    
    Beads const &
    get_beads(
        BasepairOP const &);
    
    Beads const &
    get_beads() {
        ResidueOPs res;
        beads_ = structure_->get_beads(res);
        return beads_;
    }

    inline
    void
    transform(
        Transform const & t) {
        
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
    move(
        Point const & p) { structure_->move(p); }

    String const
    to_str();
    
    String const
    to_pdb_str();
    
    void
    to_pdb(
        String const);
    
    MotifStateOP
    get_state();
    
public: //wrappers from structure
    
    inline
    AtomOPs const
    atoms() { return structure_->atoms(); }
    
    inline
    ResidueOPs const
    residues() { return structure_->residues(); }
    
    ChainOPs const &
    chains() { return structure_->chains(); }
    
    inline
    ResidueOP const
    get_residue(
        int num,
        String const & chain_id,
        String const & i_code) {
        return structure_->get_residue(num, chain_id, i_code);
    }
    
    inline
    ResidueOP const
    get_residue(
        Uuid const & uuid) {
        return structure_->get_residue(uuid);
    }
    
    void
    new_res_uuids();
    
public: //wrappers from secondary structure
    
    String
    sequence() { return secondary_structure_->sequence(); }
    
    String
    dot_bracket() { return secondary_structure_->dot_bracket(); }
    
public: //getters
    inline
    BasepairOPs const &
    ends() const { return ends_; }
    
    inline
    Strings const &
    end_ids() const { return end_ids_; }
    
    inline
    String const &
    name() const { return name_; }
    
    inline
    String const &
    path() const { return path_; }
    
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
    
    inline
    sstruct::SecondaryStructureOP const &
    secondary_structure() { return secondary_structure_; }
    
    inline
    int const &
    block_end_add() { return block_end_add_; }
    
    
public: // setters
    
    inline
    void
    mtype(MotifType const & mtype) { mtype_ = mtype; }
    
    inline
    void
    name(String const & nname) { name_ = nname; }
    
    inline
    void
    path(String const & npath) { path_ = npath; }
    
    inline
    void
    score(float const & nscore) { score_ = nscore; }
    
    inline
    void
    end_ids(Strings const & end_ids) {
        end_ids_ = end_ids;
    }
    
    inline
    void
    secondary_structure(
        sstruct::SecondaryStructureOP const & ss) {
        secondary_structure_ = ss;
    }
    
    inline
    void
    structure(
        StructureOP const & s) {
        structure_ = s;
    }
    
    inline
    void
    ends(
        BasepairOPs const & ends) { ends_ = ends; }
    
protected:
    Beads beads_;
    float score_;
    MotifType mtype_;
    BasepairOPs basepairs_, ends_;
    String path_, name_;
    Strings end_ids_;
    StructureOP structure_;
    sstruct::SecondaryStructureOP secondary_structure_;
    int block_end_add_;
};


typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP>   MotifOPs;

void
align_motif(
    BasepairStateOP const &,
    BasepairOP const &,
    MotifOP &);

MotifOP
get_aligned_motif(
    BasepairOP const &,
    BasepairOP const &,
    MotifOP const &);


Motif
ref_motif();

MotifOP
file_to_motif(
    String const &);







#endif /* defined(__RNAMake__motif__) */
