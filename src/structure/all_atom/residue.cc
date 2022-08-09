//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>

// RNAMake Headers
#include <base/exception.hpp>
#include <structure/all_atom/residue.h>
#include <util/uuid.h>

namespace structure::all_atom {

/// @brief - finds the center of the residue
math::Vector3 center_of_atoms(const Atoms &atoms) {
  auto center = math::Vector3();
  for (auto const &a : atoms) {
    center += a.get_coords();
  }
  return center / float(atoms.size());
}

/// @brief - returns the center of the residue
math::Vector3 Residue::get_center() const { return center_of_atoms(_atoms); }

/// @brief - returns the x-coordinate of the center
double Residue::get_center_x() const { return center_of_atoms(_atoms).get_x(); }

/// @brief - returns the y-coordinate of the center
double Residue::get_center_y() const { return center_of_atoms(_atoms).get_y(); }

/// @brief - returns the z-coordinate of the center
double Residue::get_center_z() const { return center_of_atoms(_atoms).get_z(); }

/// @brief - builds RNA
void Residue::_build_beads() {
  if (_rtype == structure::base::ResidueType::RNA) {
    _build_beads_RNA();
  } else if (_rtype == structure::base::ResidueType::PROTEIN) {
    _beads.push_back(util::Bead(get_coords("CA"), util::BeadType::CALPHA));
  } else {
    _beads.push_back(util::Bead(get_center(), util::BeadType::MCENTER));
  }
}

/// @brief - builds a strand of RNA
void Residue::_build_beads_RNA() {
  std::vector<Atom const *> phos_atoms, sugar_atoms, base_atoms;
  int i = -1;
  for (auto const &a : _atoms) {
    i++;
    if (a.get_name() == "P" || a.get_name() == "OP1" || a.get_name() == "OP2") {
      phos_atoms.push_back(&a);
    } else if (a.get_name().find('\'') != -1) {
      sugar_atoms.push_back(&a);
    } else {
      base_atoms.push_back(&a);
    }
  }
  auto get_center = [&](std::vector<Atom const *> const &atom_ptrs) {
    auto center = math::Vector3();
    for (auto a : atom_ptrs) {
      center = center + a->get_coords();
    }
    center = center / float(atom_ptrs.size());
    return center;
  };
  if (!phos_atoms.empty()) {
    _beads.push_back(util::Bead(get_center(phos_atoms), util::BeadType::PHOS));
  }
  if (!sugar_atoms.empty()) {
    _beads.push_back(
        util::Bead(get_center(sugar_atoms), util::BeadType::SUGAR));
  }
  if (!base_atoms.empty()) {
    _beads.push_back(util::Bead(get_center(base_atoms), util::BeadType::BASE));
  }
}

/// @brief - retrieves a residue from a .dat file (or other line/string of text
/// in a file)
Residue get_residue_from_str(const String &s) {
  Strings spl = ::base::string::split(s, ",");
  if (spl.size() < 6) {
    String msg = "String is too short! Not enough arguments";
    ::base::log_and_throw<::base::InputException>(msg);
  }
  // the name of the residue is the second space in the .dat file
  char name = spl[1][0];
  // num is the third space in the .dat file
  int num = std::stoi(spl[2]);
  // chain_id is the 4th space in the .dat file
  String chain_id = spl[3];
  // i_code will always just be a space
  char i_code = ' ';
  // uuid will be a randomly-generated thing
  util::Uuid uuid = util::generate_uuid();
  Atoms atoms;
  int i = 5;
  while (i < spl.size()) {
    if (spl[i].size() <= 1) {
      i++;
      continue;
    }
    atoms.push_back(Atom(spl[i]));
    i++;
  }
  return {name,
          num,
          chain_id,
          i_code,
          atoms,
          uuid,
          structure::base::ResidueType::RNA};
}

} // namespace structure::all_atom
