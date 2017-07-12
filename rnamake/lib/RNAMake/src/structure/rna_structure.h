//
//  rna_structure.h
//  RNAMake
//
//  Created by Joseph Yesselman on 12/2/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef RNAMAKE_STRUCTURE_RNA_STRUCTURE_H_
#define RNAMAKE_STRUCTURE_RNA_STRUCTURE_H_

#include <stdio.h>
#include <memory>
#include <algorithm>

//RNAMake
#include "primitives/rna_structure.h"
#include "structure/structure.h"
#include "structure/basepair.h"
#include "structure/residue.h"

class RNAStructure : public primitives::RNAStructure<Basepair, Structure, Chain, Residue> {
public:
    RNAStructure(
            StructureOP const &,
            BasepairOPs const &,
            BasepairOPs const &,
            SimpleStringOPs const &,
            SimpleStringOP const &,
            int,
            SimpleStringOP const &);

    RNAStructure(
            StructureOP const &,
            BasepairOPs const &,
            BasepairOPs const &,
            SimpleStringOPs const &,
            SimpleStringOP const &,
            int,
            SimpleStringOP const &,
            Beads const &);

    RNAStructure(
            RNAStructure const & rs,
            int new_uuid=0);

    RNAStructure(
            String const &,
            ResidueTypeSet const &);

    virtual
    ~RNAStructure() {}

public:
    String
    to_str();

    /*
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

*/

public: //getters
    inline
    String
    dot_bracket() { return dot_bracket_->to_str(); }

protected:
    SimpleStringOP dot_bracket_;
    Beads protein_beads_;
    int block_end_add_;

};


typedef std::shared_ptr<RNAStructure> RNAStructureOP;


/*
std::shared_ptr<BasepairOPs>
end_from_basepairs(
    StructureOP const &,
    BasepairOPs const &);

std::shared_ptr<BasepairOPs>
subselect_basepairs_with_res(
    ResidueOPs const &,
    BasepairOPs const &);
 */


//wrappers from primitive functions
inline
BasepairOPs
ends_from_basepairs(
        Structure const & s,
        BasepairOPs const & bps) {
    return primitives::ends_from_basepairs<Basepair, Structure>(s, bps);
}

inline
String
assign_end_id(
        Structure const & s,
        BasepairOPs const & bps,
        BasepairOPs const & ends,
        BasepairOP const & end) {
    return primitives::assign_end_id<Basepair, Structure, Chain, Residue>(s, bps, ends, end); }

RNAStructureOP
rna_structure_from_pdb(
        String const &,
        ResidueTypeSet const &);

BasepairOPs
basepairs_from_x3dna(
        String const &,
        Structure const &);

bool
are_rna_strucs_equal(
        RNAStructure const & rs1,
        RNAStructure const & rs2,
        bool check_uuid = 1);

#endif /* defined(RNAMAKE_STRUCTURE_RNA_STRUCTURE_H_) */


















