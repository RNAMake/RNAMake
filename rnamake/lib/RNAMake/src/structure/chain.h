//
//  chain.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__chain__
#define __RNAMake__chain__

#include <stdio.h>
#include <algorithm>

//RNAMake Headers
#include "base/types.h"
#include "structure/chain.fwd.h"
#include "structure/residue.h"
#include "structure/residue_type_set.h"

class Chain {
public:
    Chain() {}
    Chain(
        ResidueOPs const & residues):
        residues_ ( residues )
    {}
    
    Chain
    copy() const;
    
    ~Chain() {}

public:
    
    inline
    ResidueOP const &
    first() { return residues_[0]; }
    
    inline
    ResidueOP const &
    last() { return residues_.back(); }
    
    inline
    ChainOP
    subchain(int start, int end) {
        if(start < 0) { throw "start cannot be less then 0"; }
        if(start == end) { throw "start and end cannot be the same value"; }
        if(end > residues().size()) { throw "end is greater then chain length"; }
        return ChainOP(new Chain(ResidueOPs(residues_.begin() + start, residues_.begin() + end)));
    }
    
    inline
    ChainOP
    subchain(
        ResidueOP const & r1,
        ResidueOP const & r2) {
        int start = (int)(std::find(residues_.begin(), residues_.end(), r1) - residues_.begin());
        int end   = (int)(std::find(residues_.begin(), residues_.end(), r2) - residues_.begin());
        if( start > end) {
            int temp = start;
            start = end;
            end = temp;
        }
        return subchain(start, end);
    }
    
    String
    to_str() const;
    
    String
    to_pdb_str(int &) const;
    
    void
    to_pdb(String const) const;
    
public: //getters
    
    inline
    ResidueOPs &
    residues() { return residues_; }
    
private:
    ResidueOPs residues_;
};

Chain
str_to_chain(
    String const &,
    ResidueTypeSet const & );



#endif /* defined(__RNAMake__chain__) */
