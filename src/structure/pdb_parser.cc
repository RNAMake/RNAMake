//
//  pdb_parser.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 5/12/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <map>

// RNAMake Headers
//#include <base/file_io.h>
#include <base/log.hpp>
#include <base/string.hpp>
#include <math/vector_3.hpp>
#include <structure/pdb_parser.h>
#include <structure/residue.h>

namespace structure {

void PDBParser::_parse_atoms_from_pdb_file(
    String const &pdb_file, std::map<String, Atoms> &atoms) const {

  String startswith;
  String atom_name, res_name, res_num, chain_id, i_code;
  String sx, sy, sz;
  math::Vector3 atom_coords;

  auto lines = base::get_lines_from_file(pdb_file);
  auto key = String("");
  auto found = false;
  for (auto const &line : lines) {
    if (line.size() < 6) {
      continue;
    }
    startswith = line.substr(0, 6);
    if (startswith == ("ATOM  ") || startswith == ("HETATM")) {
      atom_name = line.substr(12, 4);
      atom_name = base::string::trim(atom_name);

      // do not save hydrogen atoms
      if (atom_name[0] == 'H') {
        continue;
      }

      if (_atom_name_corrections.find(atom_name) !=
          _atom_name_corrections.end()) {
        atom_name = _atom_name_corrections.at(atom_name);
      }

      res_name = line.substr(17, 4);
      res_name = base::string::trim(res_name);
      // do not save water
      if (res_name == "HOH") {
        continue;
      }

      // do not save ions
      if (_ions.find(res_name) != _ions.end()) {
        continue;
      }

      chain_id = line.substr(21, 1);
      i_code = line.substr(26, 1);

      sx = line.substr(30, 8);
      sy = line.substr(38, 8);
      sz = line.substr(46, 8);
      sx = base::string::trim(sx);
      sy = base::string::trim(sy);
      sz = base::string::trim(sz);

      atom_coords = math::Vector3(std::stod(sx), std::stod(sy), std::stod(sz));

      res_num = line.substr(22, 4);
      res_num = base::string::trim(res_num);

      key = res_name + "|" + res_num + "|" + chain_id + "|" + i_code;
      found = false;
      if (atoms.find(key) == atoms.end()) {
        atoms[key] = Atoms();
      }
      for (auto const &a : atoms[key]) {
        if (a.get_name() == atom_name) {
          found = true;
          break;
        }
      }
      if (found) {
        continue;
      }
      atoms[key].push_back(Atom(atom_name, atom_coords));

    }

    else if (startswith == "ENDMDL" || startswith.substr(0, 3) == "END") {
      break;
    }
  }
}

ResidueOP PDBParser::_setup_ref_residue(String const &pdb_file) {
  auto atoms = std::map<String, Atoms>();
  _parse_atoms_from_pdb_file(pdb_file, atoms);
  auto key = atoms.begin()->first;
  auto &res_atoms = atoms.begin()->second;
  auto spl = base::string::split(key, "|");
  auto res_name = spl[0][0];
  auto res_type = _rts.get_residue_type(spl[0]);
  auto res_num = std::stoi(spl[1]);
  auto chain_id = spl[2][0];
  auto i_code = spl[3][0];

  return std::make_shared<Residue>(res_name, res_num, chain_id, i_code,
                                   res_type, res_atoms, util::Uuid());
}

math::Matrix3x3 PDBParser::_get_res_ref_frame(ResidueCOP r) const {
  auto vec1 = math::Vector3();
  auto vec2 = math::Vector3();
  if (r->get_name() == 'A' || r->get_name() == 'G') {
    vec1 = (r->get_coords("N9") - r->get_coords("C1'")).normalize();
    vec2 =
        (r->get_coords("N9") - r->get_bead(util::BeadType::BASE).get_center())
            .normalize();
  } else {
    vec1 = (r->get_coords("N1") - r->get_coords("C1'")).normalize();
    vec2 =
        (r->get_coords("N1") - r->get_bead(util::BeadType::BASE).get_center())
            .normalize();
  }
  auto cross = vec1.cross(vec2);
  auto m = math::Matrix3x3(vec1.get_x(), vec1.get_y(), vec1.get_z(), vec2.get_x(),
                        vec2.get_y(), vec2.get_z(), cross.get_x(),
                        cross.get_y(), cross.get_z());
  m.unitarize();
  return m;
}

math::Matrix3x3
PDBParser::_get_res_ref_frame_from_atoms(std::vector<Atom const *> const &atoms,
                                         ResidueTypeCOP res_type) const {

  auto vec1 = math::Vector3();
  auto vec2 = math::Vector3();
  auto base_center = math::Vector3();
  int count = 0;
  for (int i = 12; i < res_type->get_num_atoms(); i++) {
    if (atoms[i] == nullptr) {
      continue;
    }
    base_center += atoms[i]->get_coords();
    count += 1;
  }
  base_center /= float(count);
  auto c1p_atom = atoms[res_type->get_atom_index("C1'")];

  if (res_type->get_short_name() == 'A' || res_type->get_short_name() == 'G') {
    auto n9_atom = atoms[res_type->get_atom_index("N9")];
    vec1 = (n9_atom->get_coords() - c1p_atom->get_coords()).normalize();
    vec2 = (n9_atom->get_coords() - base_center).normalize();
  } else {
    auto n1_atom = atoms[res_type->get_atom_index("N1")];
    vec1 = (n1_atom->get_coords() - c1p_atom->get_coords()).normalize();
    vec2 = (n1_atom->get_coords() - base_center).normalize();
  }

  auto cross = vec1.cross(vec2);
  auto m = math::Matrix3x3(vec1.get_x(), vec1.get_y(), vec1.get_z(), vec2.get_x(),
                        vec2.get_y(), vec2.get_z(), cross.get_x(),
                        cross.get_y(), cross.get_z());
  m.unitarize();
  return m;
}

bool PDBParser::_replace_missing_phosphate_backbone(
    std::vector<Atom const *> &atoms, ResidueTypeCOP res_type) const {

  auto ref_res = _ref_residues.at(res_type->get_name());

  // if these atoms do not exist cannot build res ref frame
  if (atoms[res_type->get_atom_index("C1'")] == nullptr) {
    return false;
  }
  if (res_type->get_short_name() == 'A' || res_type->get_short_name() == 'G') {
    if (atoms[res_type->get_atom_index("N9")] == nullptr) {
      return false;
    }
  } else {
    if (atoms[res_type->get_atom_index("N1")] == nullptr) {
      return false;
    }
  }

  auto ref_frame_1 = _get_res_ref_frame_from_atoms(atoms, res_type);
  auto ref_frame_2 = _get_res_ref_frame(ref_res);
  auto rot = dot(ref_frame_1.transpose(), ref_frame_2);
  auto r_t = rot.transpose();
  auto t = -ref_res->get_center();
  auto c4p_atom = atoms[res_type->get_atom_index("C4'")];
  if (c4p_atom == nullptr) {
    return false;
  }

  ref_res->transform(r_t, t);
  ref_res->move(c4p_atom->get_coords() - ref_res->get_coords("C4'"));

  for (int i = 0; i < 5; i++) {
    atoms[i] = new Atom(ref_res->get_atom(i).get_name(),
                        ref_res->get_atom(i).get_coords());
  }
  return true;
}

ResidueOP PDBParser::_setup_residue(String const &key, Atoms const &atoms,
                                    ResidueTypeCOP res_type) const {

  auto spl = base::string::split(key, "|");
  auto atom_ptrs = std::vector<Atom const *>(res_type->get_num_atoms());

  for (auto const &a : atoms) {
    // not a valid atom for this residue
    if (!res_type->is_valid_atom_name(a.get_name())) {
      LOGW << a.get_name() + " does not belong to residue " +
                  res_type->get_name() + ": IGNORING!";
      continue;
    }
    auto index = res_type->get_atom_index(a.get_name());
    atom_ptrs[index] = &a;
  }

  if (res_type->get_set_type() == SetType::RNA) {
    auto i = 1;
    auto missing_phosphate = false;
    for (auto const &a : atom_ptrs) {
      i++;
      if (i < 5 && a == nullptr) {
        missing_phosphate = 1;
      }
    }

    if (missing_phosphate) {
      if (!_replace_missing_phosphate_backbone(atom_ptrs, res_type)) {
        LOGW << "tried to fill in missing phosphate backbone for residue " +
                    spl[0] + " " + spl[1];
      }
    }
  }

  auto reordered_atoms = Atoms();
  auto i = 0;
  for (auto const &a_ptr : atom_ptrs) {
    if (a_ptr == nullptr) {
      LOGW << "cannot setup residue: " + spl[1] + " " + spl[0] +
                  " missing atom: " + res_type->get_atom_name_at_pos(i) +
                  " skipping this residue!";
      return nullptr;
    }
    reordered_atoms.push_back(std::move(*a_ptr));
    i++;
  }

  auto res_name = spl[0][0];
  auto res_num = std::stoi(spl[1]);
  auto chain_id = spl[2][0];
  auto i_code = spl[3][0];

  return std::make_shared<Residue>(res_name, res_num, chain_id, i_code,
                                   res_type, reordered_atoms, util::Uuid());
}

PDBParserResiduesOP PDBParser::parse(String const &pdb_file) const {

  auto residues = std::make_shared<PDBParserResidues>();
  auto atoms = std::map<String, Atoms>();
  _parse_atoms_from_pdb_file(pdb_file, atoms);
  for (auto const &kv : atoms) {
    auto spl = base::string::split(kv.first, "|");
    auto has_res_type = _rts.contains_residue_type(spl[0]);
    if (has_res_type) {
      auto res_type = _rts.get_residue_type(spl[0]);
      if (res_type->get_set_type() == SetType::RNA) {
        auto r = _setup_residue(kv.first, kv.second, res_type);
        if (r == nullptr) {
          continue;
        }
        residues->RNA_residues.push_back(r);
      } else if (res_type->get_set_type() == SetType::PROTEIN) {
        auto r = _setup_residue(kv.first, kv.second, res_type);
        if (r == nullptr) {
          continue;
        }
        residues->protein_residues.push_back(r);
      }
    } else {
      auto atom_names = Strings();
      for (auto const &a : kv.second) {
        atom_names.push_back(a.get_name());
      }
      auto res_type = get_new_residue_type(spl[0], atom_names);
      residues->small_molecule_residues.push_back(
          _setup_residue(kv.first, kv.second, res_type));
    }
  }
  return residues;
}
} // namespace structure
