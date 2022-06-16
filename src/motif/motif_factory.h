//
//  motif_factory.h
//  RNAMake
//
//  Created by Joseph Yesselman on 8/1/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__motif_factory__
#define __RNAMake__motif_factory__

#include <filesystem>

#include <stdio.h>
#include <util/x3dna.h>

// RNAMake Headers
#include "base/env_manager.h"
#include "base/file_io.h"
#include "base/settings.h"
#include "base/types.hpp"
#include "motif/motif.h"
#include "motif/motif_scorer.h"
#include "motif/motif_to_secondary_structure.h"
#include "structure/pdb_parser.h"

namespace motif {

/*
 * Exception for Motif Factory
 */
class MotifFactoryException : public std::runtime_error {
public:
  /**
   * Standard constructor for MotifFactoryException
   * @param   message   Error message for Motif Factory
   */
  MotifFactoryException(String const &message) : std::runtime_error(message) {}
};

class MotifFactory {
public:
  MotifFactory()
      : parser_(MotiftoSecondaryStructure()),
        pdb_parser_(structure::PDBParser()) {
    // make sure x3dna is set
    util::X3dna::set_envs();
    auto path = base::motif_dirs() + "ref.motif";
    ref_motif_ = file_to_motif(path);
    path = base::motif_dirs() + "base.motif";
    base_motif_ = file_to_motif(path);
    base_motif_->get_beads(base_motif_->ends()[0]);
    added_helix_ = std::make_shared<Motif>(*base_motif_);
    clash_radius_ = 2.9;
  }

  ~MotifFactory() {}

public:
  MotifOP motif_from_file(String const &path, bool rebuild_x3dna = true,
                          bool include_protein = false,
                          int force_num_chains = -1);

  MotifOP get_oriented_motif(MotifOP const &, int);

  MotifOP motif_from_res(structure::ResidueOPs &,
                         structure::BasepairOPs const &);

  MotifOP motif_from_bps(structure::BasepairOPs const &);

  MotifOP can_align_motif_to_end(MotifOP const &, int);

  MotifOP align_motif_to_common_frame(MotifOP const &, int);

  void standardize_rna_structure_ends(MotifOP &);

  structure::BasepairOPs _setup_basepair_ends(structure::StructureOP const &,
                                              structure::BasepairOPs const &);

  void _setup_secondary_structure(MotifOP &);

public:
  inline MotifOP const &added_helix() { return added_helix_; }

  inline MotifOP const &ref_motif() { return ref_motif_; }

private:
  void _standardize_motif(MotifOP &);

  structure::BasepairOPs _setup_basepairs(String const &,
                                          structure::StructureOP const &, bool);

  void _align_chains(MotifOP &);

  void _align_ends(MotifOP &);

  int _steric_clash(MotifOP const &, MotifOP const &);

  int _steric_clash_count(MotifOP const &, MotifOP const &);

  int _bead_overlap(MotifOP const &, MotifOP const &);

private: // new functions to add to new code
  structure::StructureOP
  _get_reduced_chain_num_structure(structure::Structure const &, int);

private:
  structure::PDBParser pdb_parser_;
  MotiftoSecondaryStructure parser_;
  MotifScorer scorer_;
  MotifOP ref_motif_, base_motif_, added_helix_;
  float clash_radius_;
};

} // namespace motif

#endif /* defined(__RNAMake__motif_factory__) */
