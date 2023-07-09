//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <fstream>

//RNAMake Headers
#include "structure/residue.h"

namespace structure {

math::Point
center( AtomOPs const & atoms) {
    assert(atoms.size() > 0);
    auto center = math::Point(0, 0, 0);
    for (auto const & a : atoms) {
        if (a == nullptr) { continue; }
        center += a->coords();
    }
    return center / float(atoms.size());
}


void
Residue::setup_atoms(
        AtomOPs & atoms) {
    _atoms = AtomOPs(_rtype.size());
    int count = 0;
    for (auto & a : atoms) {
        if (a == nullptr) { continue; }
        auto name_change = _rtype.get_correct_atom_name(*a);
        //check for misnamed atoms
        if (name_change.length() != 0) { a->name(name_change); }
        auto pos = _rtype.atom_pos_by_name(a->name());
        if (pos == -1) { continue; }
        _atoms[pos] = a;
        count++;
    }

    if (count == 0) { throw ResidueException("Residue has zero atoms something wrong"); }
}

structure::Beads
Residue::get_beads() const {
    AtomOPs phos_atoms, sugar_atoms, base_atoms;
    int i = -1;
    for (auto const & a : _atoms) {
        i++;
        if (a == nullptr) { continue; }
        if (i < 3) { phos_atoms.push_back(a); }
        else if (i < 12) { sugar_atoms.push_back(a); }
        else { base_atoms.push_back(a); }
    }
    auto beads = Beads();
    if (phos_atoms.size() > 0) {
        beads.push_back(Bead(::structure::center(phos_atoms), BeadType::PHOS));
    }
    if (sugar_atoms.size() > 0) {
        beads.push_back(Bead(::structure::center(sugar_atoms), BeadType::SUGAR));
    }
    if (base_atoms.size() > 0) {
        beads.push_back(Bead(::structure::center(base_atoms), BeadType::BASE));
    }
    return beads;
}


String
Residue::to_str() const {
    auto s = String();
    s = _rtype.name() + "," + _name + "," + std::to_string(_num) + "," + _chain_id + "," + _i_code + ",";
    for (auto const & a : _atoms) {
        if (a == nullptr) { s += "N,"; }
        else { s += a->to_str() + ","; }
    }
    return s;
}

String
Residue::to_pdb_str(
        int & acount,
        int rnum,
        String const & chain_id) const {

    auto num = _num;
    if (rnum != -1) {
        num = rnum;
    }
    auto cid = _chain_id;
    if (chain_id != "") {
        cid = chain_id;
    }

    auto s = String();
    for (auto const & a : _atoms) {
        if (a == nullptr) { continue; }
        char buffer[200];
        std::sprintf(buffer, "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", "ATOM", acount,
                     a->name().c_str(), "", short_name().c_str(), cid.c_str(), num, "", a->coords()[0], a->coords()[1],
                     a->coords()[2], 1.00, 0.00, "", "");
        s += String(buffer);
        acount++;
    }
    return s;
}

void
Residue::to_pdb(
        String const fname) {
    std::ofstream out;
    out.open(fname.c_str());
    int i = 1;
    auto s = to_pdb_str(i);
    out << s << std::endl;
    out.close();
}

}



