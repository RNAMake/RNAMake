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
#include <base/settings.h>
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
              _rts(rts),
              _atom_name_corrections(std::map<String, String>()),
              _ions(std::map<String, int>()) {

          _ions["MG"] = 1;
          _ions["K"] = 1;
          _ions["BR"] = 1;
          _ions["ZN"] = 1;
          _ions["MN"] = 1;
          _ions["NA"] = 1;
          _ions["CL"] = 1;
          _ions["PO4"] = 1;
          _ions["SO4"] = 1;
          _ions["NH4"] = 1;
          _ions["NO3"] = 1;

          _atom_name_corrections["O1P"] = "OP1";
          _atom_name_corrections["O2P"] = "OP2";

          auto path = base::resources_path() + "/ideal_residues/";
          auto names = Strings{"ADE", "CYT", "GUA", "URA"};
          for(auto const & name : names) {
              _ref_residues[name] = _setup_ref_residue(path + name + ".pdb");
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
      ResidueTypeSet const & _rts;
      std::map<String, ResidueOP> _ref_residues;
      std::map<String, String> _atom_name_corrections;
      std::map<String, int> _ions;


  };

}
#endif /* defined(__RNAMake__structure_pdb_parser__) */