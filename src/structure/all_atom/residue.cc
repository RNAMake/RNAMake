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


String Residue::get_str() const {
  std::stringstream ss;
  // ss << _res_type->get_name() << "," << _name << "," << _num << "," <<
  //     _chain_id
  //    << "," << _i_code << ",";
  for (auto const &a : _atoms) {
  //  ss << a.get_str() + ",";
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
