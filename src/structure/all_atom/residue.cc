//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>

// RNAMake Headers
#include <structure/all_atom/residue.h>

namespace structure::all_atom {

math::Vector3 center_of_atoms(const Atoms &atoms) {
  auto center = math::Vector3();
  for (auto const &a : atoms) {
    center += a.get_coords();
  }
  return center / float(atoms.size());
}

void Residue::_build_beads() {
  if (_rtype == ResidueType::RNA) {
    _build_beads_RNA();
  } else if (_rtype == ResidueType::PROTEIN) {
    _beads.push_back(util::Bead(get_coords("CA"), util::BeadType::CALPHA));
  } else {
    _beads.push_back(util::Bead(get_center(), util::BeadType::MCENTER));
  }
}

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

String Residue::get_str() const {
  std::stringstream ss;
  // ss << _res_type->get_name() << "," << _name << "," << _num << "," <<
  //     _chain_id
  //    << "," << _i_code << ",";
  for (auto const &a : _atoms) {
    ss << a.get_str() + ",";
  }
  return ss.str();
}

String Residue::get_bead_pdb_str(int &acount, int rnum, char chain_id) const {
  auto s = String();
  auto bead_name = String("C");
  for (auto const &b : _beads) {
    char buffer[200];
    auto c = b.get_center();
    if (b.get_type() == util::BeadType::BASE) {
      bead_name = "C";
    } else if (b.get_type() == util::BeadType::SUGAR) {
      bead_name = "N";
    } else if (b.get_type() == util::BeadType::PHOS) {
      bead_name = "P";
    }
    /*
    std::sprintf(buffer,
                 "%-6s%5d %-4s%1s%-4c%1c%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f     "
                 " %4s%2s\n",
                 "ATOM", acount, bead_name.c_str(), "", _name, chain_id, rnum,
                 "", c.get_x(), c.get_y(), c.get_z(), 1.00, 0.00, "", "");
    s += String(buffer); */
    acount++;
  }
  return s;
}


/*
void Residue::write_pdb(String const fname) {
  std::ofstream out;
  out.open(fname.c_str());
  out << get_pdb_str(1) << std::endl;
  out.close();
} */

} // namespace structure::all_atom
