//
//  rna_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__sec_rna_structure__
#define __RNAMake__sec_rna_structure__

#include <stdio.h>
#include <algorithm>


//RNAMake
#include "secondary_structure/basepair.h"
#include "secondary_structure/structure.h"

namespace secondary_structure {


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
    
    RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends,
        Strings const & end_ids,
        String const & name,
        String const & path,
        float score):
    structure_(structure),
    basepairs_(basepairs),
    ends_(ends),
    name_(name),
    path_(path),
    score_(score),
    end_ids_(end_ids)
    {}
  
    inline
    RNAStructure(
        RNAStructure const & rs):
    structure_( std::make_shared<Structure>(*rs.structure_)),
    basepairs_( BasepairOPs(rs.basepairs_.size()) ),
    ends_( BasepairOPs(rs.ends_.size()) )
    {
    
        int i = 0;
        for(auto const & bp : rs.basepairs_) {
            auto new_bp = std::make_shared<Basepair>(structure_->get_residue(bp->res1()->uuid()),
                                                     structure_->get_residue(bp->res2()->uuid()),
                                                     bp->uuid());
            basepairs_[i] = new_bp;
            i++;
        }
        
        i = 0;
        for(auto const & end : rs.ends_) {
            int pos = (int)(std::find(rs.basepairs_.begin(),
                                      rs.basepairs_.end(), end) - rs.basepairs_.begin());
            ends_[i] = basepairs_[pos];
            i++;
            
        }
        
        name_    = rs.name_;
        path_    = rs.path_;
        end_ids_ = rs.end_ids_;
    
    }
    
    virtual
    ~RNAStructure() {}
    
public:
    
    BasepairOPs
    get_basepair(
        String const &);
    
    BasepairOPs
    get_basepair(
        util::Uuid const &);
    
    BasepairOPs
    get_basepair(
        ResidueOP const &,
        ResidueOP const &);
    
    BasepairOPs
    get_basepair(
        util::Uuid const &,
        util::Uuid const &);

public:
    BasepairOP
    get_end(
            String const &);

    virtual
    void
    replace_sequence(
        String const &);
    
public: //wrappers for structure
    
    inline
    ResidueOP
    get_residue(
        int num,
        String const & chain_id,
        String const & i_code) {
        return structure_->get_residue(num, chain_id, i_code);
    }
    
    inline
    ResidueOP
    get_residue(
        util::Uuid const & uuid) {
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
    
    inline
    StructureOP
    structure() { return structure_; }
    
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
