//
//  motif.h
//  RNAMake
//
//  Created by Joseph Yesselman on 7/31/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__ss_motif__
#define __RNAMake__ss_motif__

#include <stdio.h>

#include "util/motif_type.h"
#include "secondary_structure/chain.h"
#include "secondary_structure/basepair.h"
#include "secondary_structure/rna_structure.h"

namespace sstruct {
 
class Motif : public RNAStructure {
public:
    
    Motif():
    RNAStructure()
    {}
    
    Motif(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    RNAStructure(structure, basepairs, ends)
    {}


    Motif(
        Motif const & motif) {
        structure_ = std::make_shared<Structure>(*motif.structure_);
    }
    
    
    ~Motif() {}
    
 
public: //getters
    
    inline
    MotifType const &
    mtype() { return mtype_; }
    

public: //setters
    
    inline
    void
    mtype(MotifType const & mtype) { mtype_ = mtype; }
    
    
private:
    MotifType mtype_;
    
};
    
typedef std::shared_ptr<Motif> MotifOP;
typedef std::vector<MotifOP>   MotifOPs;
    
}

#endif /* defined(__RNAMake__motif__) */
