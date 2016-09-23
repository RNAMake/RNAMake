//
//  motif_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_factory__
#define __RNAMake__motif_factory__

#include <stdio.h>

//RNAMake Headers
#include "base/types.h"
#include "base/settings.h"
#include "base/file_io.h"
#include "structure/pdb_parser.h"
#include "motif/motif_to_secondary_structure.h"
#include "motif/motif_scorer.h"
#include "motif/motif.h"


/*
 * Exception for Motif Factory
 */
class MotifFactoryException : public std::runtime_error {
public:
    /**
     * Standard constructor for MotifFactoryException
     * @param   message   Error message for Motif Factory
     */
    MotifFactoryException(String const & message):
    std::runtime_error(message)
    {}
};

class MotifFactory {
public:
    MotifFactory():
    parser_(MotiftoSecondaryStructure()),
    pdb_parser_(PDBParser()) {
        auto path = motif_dirs() + "ref.motif";
        ref_motif_ = file_to_motif(path);
        path = motif_dirs() + "base.motif";
        base_motif_ = file_to_motif(path);
        base_motif_->get_beads(base_motif_->ends()[0]);
        added_helix_ = std::make_shared<Motif>(*base_motif_);
        clash_radius_ = 2.9;
    }
    
    ~MotifFactory() {}
    
public:
    
    MotifOP
    motif_from_file(
        String const & path,
        bool rebuild_x3dna = true,
        bool include_protein = false);
    
    MotifOP
    motif_from_res(
        ResidueOPs &,
        BasepairOPs const &);
    
    MotifOP
    motif_from_bps(
        BasepairOPs const &);
    
    void
    standardize_motif(
        MotifOP &);
    
    MotifOP
    can_align_motif_to_end(
        MotifOP const &,
        int);
    
    MotifOP
    align_motif_to_common_frame(
        MotifOP const &,
        int);
    
    BasepairOPs
    _setup_basepair_ends(
        StructureOP const &,
        BasepairOPs const &);
    
    void
    _setup_secondary_structure(
        MotifOP &);
    
public:
    inline
    MotifOP const &
    added_helix() { return added_helix_; }
    
    
private:
    
    BasepairOPs
    _setup_basepairs(
        String const &,
        StructureOP const &,
        bool);
    
    
    void
    _align_chains(
        MotifOP &);
    
    void
    _align_ends(
        MotifOP &);
    
    int
    _steric_clash(
        MotifOP const &,
        MotifOP const &);
    
    
private:
    MotiftoSecondaryStructure parser_;
    PDBParser pdb_parser_;
    MotifScorer scorer_;
    MotifOP ref_motif_, base_motif_, added_helix_;
    float clash_radius_;
    
    
};

#endif /* defined(__RNAMake__motif_factory__) */
