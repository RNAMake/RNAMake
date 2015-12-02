//
//  rna_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__rna_structure__
#define __RNAMake__rna_structure__

#include <stdio.h>

//RNAMake
#include "secondary_structure/basepair.h"
#include "secondary_structure/structure.h"

namespace sstruct {


class RNAStructure {
public:
    RNAStructure():
    structure_(StructureOP()),
    basepairs_(BasepairOPs()),
    ends_(BasepairOPs()),
    name_(""),
    path_(""),
    score_(0),
    end_ids_(Strings())
    {}
    
    RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    structure_(structure),
    basepairs_(basepairs),
    ends_(ends),
    name_(""),
    path_(""),
    score_(0),
    end_ids_(Strings())
    {}
    
public:
    
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
    
public: //wrappers for structure
    
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

    inline
    String
    sequence() { return structure_->sequence(); }
    
    inline
    String
    dot_bracket() { return structure_->dot_bracket(); }
    
    inline
    ChainOPs const &
    chains() { return structure_->chains(); }
    
    inline
    ResidueOPs
    residues() { return structure_->residues(); }
    
public: //getters
    
    inline
    BasepairOPs const &
    basepairs() { return basepairs_; }
    
    inline
    BasepairOPs const &
    ends() { return ends_; }
    
    inline
    String const &
    name() { return name_; }
    
    inline
    Strings const &
    end_ids() { return end_ids_; }

public: //setters
    
    inline
    void
    name(String const & name) { name_ = name; }
    
    inline
    void
    path(String const & path) { path_ = path; }
    
    inline
    void
    end_ids(Strings const & end_ids) { end_ids_ = end_ids; }
    
protected:
    StructureOP structure_;
    BasepairOPs basepairs_, ends_;
    String name_, path_;
    Strings end_ids_;
    float score_;
    
};

typedef std::shared_ptr<RNAStructure> RNAStructureOP;
    
}

#endif /* defined(__RNAMake__rna_structure__) */
