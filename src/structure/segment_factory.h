//
// Created by Joseph Yesselman on 12/20/17.
//

#ifndef RNAMAKE_NEW_SEGMENT_FACTORY_H
#define RNAMAKE_NEW_SEGMENT_FACTORY_H

#include <util/x3dna.h>
#include "structure/aligner.h"
#include "structure/segment.h"
#include "structure/pdb_parser.h"
#include "structure/basepair.h"
#include "structure/structure.h"
#include "structure/pose.h"
#include "structure/residue_type_set.h"

namespace structure {

  class SegmentFactoryException : public std::runtime_error {
  public:
      SegmentFactoryException(
              String const & message):
              std::runtime_error(message) {}
  };

  class SegmentFactory {
  public:
      struct SegmentElements {
          base::SimpleStringCOP name;
          base::SimpleStringCOPs end_ids;
          Structure rna, proteins, small_molecules;
          Basepairs basepairs;
          Indexes end_indexes;

          inline
          SegmentElements(
                  base::SimpleStringCOP nname,
                  Structure const & rna_struc,
                  Structure const & protein_struc,
                  Structure const & small_molecule_struc,
                  Basepairs const & bps,
                  Indexes const & nend_indexes):
                  name(nname),
                  end_ids(base::SimpleStringCOPs()),
                  rna(std::move(rna_struc)),
                  proteins(std::move(protein_struc)),
                  small_molecules(std::move(small_molecule_struc)),
                  basepairs(std::move(bps)),
                  end_indexes(nend_indexes) {}

          inline
          Basepair &
          get_end(
                  Index end_index) {
              return basepairs[end_indexes[end_index]];
          }

          inline
          Basepairs
          get_ends() {
              auto ends = Basepairs();
              for(auto const & ei : end_indexes) {
                  ends.push_back(basepairs[ei]);
              }
              return ends;
          }

          inline
          void
          move(
                  math::Point const & p) {
              rna.move(p);
              proteins.move(p);
              small_molecules.move(p);
              for(auto & bp : basepairs) { bp.move(p); }
          }

          inline
          void
          transform(
                  math::Matrix const & r,
                  math::Point const & t) {
              auto dummy = math::Point();
              rna.transform(r, t, dummy);
              proteins.transform(r, t, dummy);
              small_molecules.transform(r, t, dummy);
              for(auto & bp : basepairs) { bp.transform(r, t, dummy); }
          }

      };

      typedef std::shared_ptr<SegmentElements> SegmentElementsOP;

  public:
      SegmentFactory(
              ResidueTypeSet const & rts):
              rts_(rts),
              x3dna_(util::X3dna()),
              pdb_parser_(PDBParser(rts)) {
          // dont need to rebuild x3dna files for ref motifs
          x3dna_.set_rebuild_files(false);
          ref_motif_  = _setup_ref_motif();
          base_helix_ = _setup_base_helix();
          added_helix_ = std::make_shared<Segment>(*base_helix_);
          x3dna_.set_rebuild_files(true);
      }

  public:
      SegmentOP
      segment_from_pdb(
              String const & pdb_path,
              util::SegmentType segment_type = util::SegmentType::SEGMENT,
              bool rebuild_x3dna_files = true) const;

      SegmentOPs
      all_segments_from_pdb(
              String const & pdb_path,
              util::SegmentType segment_type = util::SegmentType::SEGMENT,
              bool rebuild_x3dna_files = true) const;

      SegmentOP
      segment_from_components(
              String const & name,
              Structure const & rna_struc,
              Basepairs const & basepairs,
              Structure const & proteins,
              Structure const & small_molecules,
              util::SegmentType segment_type = util::SegmentType::SEGMENT) const;

      void
      align_segment_to_ref_frame(
              Segment &) const;

  private:
      PoseOP
      _setup_ref_motif();

      SegmentOP
      _setup_base_helix();

  public:

      SegmentElementsOP
      _get_segment_elements_from_pdb(
              String const &,
              bool) const;

      void
      _check_common_segment_issues(
              SegmentElements const &,
              util::SegmentType) const;

  private:

      void
      _setup_end_ids(
              Structure const &,
              Basepairs const &,
              Indexes const &,
              base::SimpleStringCOPs &) const;

      StructureOP
      _get_aligned_structure(
              Structure const &) const;

      void
      _get_aligned_end_indexes(
              Structure const &,
              Indexes &,
              Basepairs &) const;

      void
      _remove_beads_from_end_res(
              Basepairs const &,
              Indexes const &,
              Structure &) const;

      void
      _standardize_motif_elements(
              SegmentElements &,
              Index) const;

      void
      _align_motif_elements_to_frame(
              Basepair const &,
              SegmentElements &,
              Index) const;

      void
      _align_motif_elements_back_to_org_frame(
              Basepairs const &,
              Structure const &,
              SegmentElements &) const;

      void
      _align_segment_to_frame(
              Basepair const &,
              Segment &) const;

      SegmentOP
      _get_aligned_segment(
              Basepair const &,
              Segment &) const;

      int
      _num_steric_clashes(
              SegmentElements const &,
              Segment const &) const;

      int
      _are_structures_overlaid(
              Structure const &,
              Structure const &) const;



  private:
      ResidueTypeSet const & rts_;
      Aligner aligner_;
      util::X3dna x3dna_;
      PDBParser pdb_parser_;
      PoseOP ref_motif_;
      SegmentOP base_helix_, added_helix_;
      // variables for aligning
  };

  typedef std::shared_ptr<SegmentFactory> SegmentFactoryOP;


}

#endif //RNAMAKE_NEW_SEGMENT_FACTORY_H