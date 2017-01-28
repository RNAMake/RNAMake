//
//  rna_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMAKE_LIB_RNAMAKE_SRC_STRUCTURE_RNA_STRUCTURE_H_
#define RNAMAKE_LIB_RNAMAKE_SRC_STRUCTURE_RNA_STRUCTURE_H_

#include <stdio.h>
#include <memory>
#include <algorithm>

//RNAMake
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/residue.h"


/*
 * Exception for RNA Structure
 */
class RNAStructureException : public std::runtime_error {
public:
    /**
     * Standard constructor for RNAStructureException
     * @param   message   Error message for rna structure
     */
    RNAStructureException(String const & message):
    std::runtime_error(message)
    {}
};

class RNAStructure {
public:
    RNAStructure()
    {}
    
    RNAStructure(
        StructureOP const & structure,
        BasepairOPs const & basepairs,
        BasepairOPs const & ends):
    structure_(structure),
    basepairs_(basepairs),
    ends_(ends)
    {}
    
    virtual
    ~RNAStructure() {}
    
public: // get specific basepairs
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
    
    
public: // get steric beads
    Beads const &
    get_beads(
        BasepairOPs const &);
    
    Beads const &
    get_beads(
        BasepairOP const &);
    
    inline
    Beads const &
    get_beads() {
        ResidueOPs res;
        beads_ = structure_->get_beads(res);
        return beads_;
    }
    
public: //get end information
    int
    get_end_index(
        BasepairOP const &);
    
    int
    get_end_index(
        String const &);
    
public: //output functions
    
    String const
    to_pdb_str(
        int rnumber = -1);
    
    void
    to_pdb(
        String const,
        int renumber = -1);
    
public: //wrappers from structure

    inline
    ResidueOPs const
    residues() const { return structure_->residues(); }
    
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
    float const &
    score() const { return score_; }
    
    inline
    Beads const &
    protein_beads() { return protein_beads_; }

    
public: //setters
    
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
    ends(BasepairOPs const & ends) { ends_ = ends; }
    
    inline
    void
    protein_beads(Beads const & beads) { protein_beads_ = beads; }
    
    
protected:
    StructureOP structure_;
    BasepairOPs basepairs_, ends_;
    String name_, path_;
    Strings end_ids_;
    Beads beads_, protein_beads_;
    float score_;
    
};


typedef std::shared_ptr<RNAStructure> RNAStructureOP;

std::shared_ptr<BasepairOPs>
end_from_basepairs(
    StructureOP const &,
    BasepairOPs const &);

std::shared_ptr<BasepairOPs>
subselect_basepairs_with_res(
    ResidueOPs const &,
    BasepairOPs const &);



#endif /* defined(RNAMAKE_LIB_RNAMAKE_SRC_STRUCTURE_RNA_STRUCTURE_H_) */


















