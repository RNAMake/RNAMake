//
//  pdb_parser.h
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__structure_pdb_parser__
#define __RNAMake__structure_pdb_parser__

#include <stdio.h>
#include <map>

//RNAMake Headers
#include <base/paths.h>
#include <math/xyz_matrix.h>
#include <structure/residue.h>
#include <structure/residue_type_set.h>


namespace structure {
  struct PDBParserResidues {
      ResidueOPs RNA_residues, protein_residues, small_molecule_residues;

      inline
      bool
      has_RNA() {
          if(RNA_residues.size() > 0) { return true; }
          else                        { return false; }
      }

      inline
      bool
      has_protein() {
          if(protein_residues.size() > 0) { return true; }
          else                            { return false; }
      }

      inline
      bool
      has_small_molecules() {
          if(small_molecule_residues.size() > 0) { return true; }
          else                                   { return false; }
      }
  };

  typedef std::shared_ptr<PDBParserResidues> PDBParserResiduesOP;

  class PDBParser {
  public:
      PDBParser(
              ResidueTypeSet const & rts) :
              rts_(rts),
              atom_name_corrections_(std::map<String, String>()),
              ions_(std::map<String, int>()) {

          ions_["MG"] = 1;
          ions_["K"] = 1;
          ions_["BR"] = 1;
          ions_["ZN"] = 1;
          ions_["MN"] = 1;
          ions_["NA"] = 1;
          ions_["CL"] = 1;
          ions_["PO4"] = 1;
          ions_["SO4"] = 1;
          ions_["NH4"] = 1;
          ions_["NO3"] = 1;

          atom_name_corrections_["O1P"] = "OP1";
          atom_name_corrections_["O2P"] = "OP2";

          auto path = base::resources_path() + "/ideal_residues/";
          auto names = Strings{"ADE", "CYT", "GUA", "URA"};
          for(auto const & name : names) {
              ref_residues_[name] = _setup_ref_residue(path + name + ".pdb");
          }
      }

      ~PDBParser() {}

  public:

      PDBParserResiduesOP
      parse(
              String const &) const;

  private:
      void
      _parse_atoms_from_pdb_file(
              String const &,
              std::map<String, Atoms> &) const;

      ResidueOP
      _setup_ref_residue(
              String const &);

      ResidueOP
      _setup_residue(
              String const &,
              Atoms const &,
              ResidueTypeCOP) const;

      math::Matrix
      _get_res_ref_frame(
              ResidueCOP) const;

      math::Matrix
      _get_res_ref_frame_from_atoms(
              std::vector<Atom const *> const &,
              ResidueTypeCOP) const;

      bool
      _replace_missing_phosphate_backbone(
              std::vector<Atom const *> &,
              ResidueTypeCOP) const;


  private:
      ResidueTypeSet const & rts_;
      std::map<String, ResidueOP> ref_residues_;
      std::map<String, String> atom_name_corrections_;
      std::map<String, int> ions_;


  };

}
#endif /* defined(__RNAMake__structure_pdb_parser__) */