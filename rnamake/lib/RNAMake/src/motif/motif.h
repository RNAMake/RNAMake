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

class Motif : public RNAStructure {
public:
    Motif():
    RNAStructure(),
    id_(Uuid()),
    block_end_add_(0)
    {}
    
    Motif(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    RNAStructure(structure, basepairs, ends),
    id_(Uuid()),
    block_end_add_(0)
    {}
        
    Motif(String const &,
          ResidueTypeSet const &);
    
    
    Motif(
        Motif const & m);
    
    ~Motif() {}

public:

    
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
    
      
    MotifStateOP
    get_state();
    
    
    void
    new_res_uuids();
    
public: //wrappers from secondary structure
    
    String
    sequence() { return secondary_structure_->sequence(); }
    
    String
    dot_bracket() { return secondary_structure_->dot_bracket(); }
    
public: //getters
    inline
    MotifType const &
    mtype() { return mtype_; }
    
    inline
    sstruct::MotifOP const &
    secondary_structure() { return secondary_structure_; }
    
    inline
    int const &
    block_end_add() { return block_end_add_; }
    
    inline
    Uuid const &
    id() { return id_; }
    
    
public: // setters
    
    inline
    void
    id(Uuid const & nid) { id_ = nid; }
    
    inline
    void
    mtype(MotifType const & mtype) { mtype_ = mtype; }
    
    inline
    void
    secondary_structure(
        sstruct::MotifOP const & ss) {
        secondary_structure_ = ss;
    }
    
    inline
    void
    structure(
        StructureOP const & s) {
        structure_ = s;
    }
    
   
    
private:
    MotifType mtype_;
    sstruct::MotifOP secondary_structure_;
    int block_end_add_;
    Uuid id_;
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
