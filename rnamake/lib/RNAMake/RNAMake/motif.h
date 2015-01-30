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
    
    ~Motif() {}
 
    
private:
    Beads beads_;
    float score_;
    MotifType mtype_;
    Basepairs basepairs_, ends_;
    String mdir_, name_;
    Matrices cached_rotations_;
    Structure structure_;
};

#endif /* defined(__RNAMake__motif__) */
