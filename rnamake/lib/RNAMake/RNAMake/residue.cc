//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>
#include <fstream>
#include "residue.h"
#include "util.h"

void
Residue::setup_atoms(
    AtomOPs & atoms) {
    atoms_ = AtomOPs(atoms.size(), NULL);
    for(auto & a : atoms) {
        if(a == NULL) { continue; }
        String name_change = rtype_.get_correct_atom_name(*a);
        //check for misnamed atoms
        if(name_change.length() != 0) { a->name(name_change); }
        int pos = rtype_.atom_pos_by_name(a->name());
        if(pos == -1) { continue; }
        atoms_[pos] = a;
    }
}

Beads
Residue::get_beads() {
    AtomOPs phos_atoms, sugar_atoms, base_atoms;
    int i = -1;
    for( auto const & a : atoms_) {
        i++;
        if(a == NULL) { continue; }
        if     ( i < 3 ) { phos_atoms.push_back(a);  }
        else if( i < 12) { sugar_atoms.push_back(a); }
        else             { base_atoms.push_back(a);  }
    }
    Beads beads;
    if(phos_atoms.size() > 0)  { beads.push_back(Bead(center(phos_atoms), PHOS)); }
    if(sugar_atoms.size() > 0) { beads.push_back(Bead(center(sugar_atoms), SUGAR)); }
    if(base_atoms.size() > 0)  { beads.push_back(Bead(center(base_atoms), BASE)); }
    return beads;
}

Residue
Residue::copy() const {
    Residue copied_r (rtype_, name_, num_, chain_id_, i_code_);
    AtomOPs copied_atoms(atoms_.size());
    int i = -1;
    for(auto const & a : atoms_) {
        i++;
        if(a == NULL) { continue; }
        copied_atoms[i] = AtomOP( new Atom( a->copy()));
    }
    copied_r.uuid(uuid_);    
    return copied_r;
}

String
Residue::to_str() const {
    std::stringstream ss;
    ss << rtype_.name() << "," << name_ << "," << num_ << "," << chain_id_ << "," << i_code_ << ",";
    for( auto const & a : atoms_) {
        if( a == NULL) { ss << "N,"; }
        else           { ss << a->to_str() + ","; }
    }
    return ss.str();
}

String
Residue::to_pdb_str(
    int & acount) const {

    String s;
    for(auto const & a : atoms_) {
        if( a == NULL) { continue; }
        char buffer [200];
        std::sprintf(buffer, "%-6s%5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n", "ATOM", acount, a->name().c_str(), "", short_name().c_str(), chain_id_.c_str(), num_, "",  a->coords()[0], a->coords()[1], a->coords()[2], 1.00, 0.00, "", "");
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
    String s = to_pdb_str(i);
    out << s << std::endl;
    out.close();
}






