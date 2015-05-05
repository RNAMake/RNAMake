//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <sstream>
#include <fstream>

//RNAMake Headers
#include "structure/residue.h"

Point
center(
    AtomOPs const & atoms) {
    assert(atoms.size() > 0);
    Point center(0, 0, 0);
    for(auto const & a : atoms) {
        if(a == nullptr) {
            continue;
        }
        center += a->coords();
    }
    return center / float(atoms.size());
}


void
Residue::setup_atoms(
    AtomOPs & atoms) {
    atoms_ = AtomOPs(atoms.size());
    int count = 0;
    for(auto & a : atoms) {
        if(a == nullptr) { continue; }
        String name_change = rtype_.get_correct_atom_name(*a);
        //check for misnamed atoms
        if(name_change.length() != 0) { a->name(name_change); }
        int pos = rtype_.atom_pos_by_name(a->name());
        if(pos == -1) { continue; }
        atoms_[pos] = a;
        count++;
    }
    
    if(count == 0) { throw "Residue has zero atoms something wrong"; }
}

Beads
Residue::get_beads() const {
    AtomOPs phos_atoms, sugar_atoms, base_atoms;
    int i = -1;
    for( auto const & a : atoms_) {
        i++;
        if(a  == nullptr) { continue; }
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
    copied_r.atoms_ = AtomOPs(atoms_.size());
    int i = -1;
    for(auto const & a : atoms_) {
        i++;
        if(a  == nullptr) { continue; }
        copied_r.atoms_[i] = AtomOP( new Atom( a->copy()));
    }
    copied_r.uuid_ = uuid_;
    return copied_r;
}

String
Residue::to_str() const {
    std::stringstream ss;
    ss << rtype_.name() << "," << name_ << "," << num_ << "," << chain_id_ << "," << i_code_ << ",";
    for( auto const & a : atoms_) {
        if( a == nullptr) { ss << "N,"; }
        else           { ss << a->to_str() + ","; }
    }
    return ss.str();
}

String
Residue::to_pdb_str(
    int & acount) const {

    String s;
    for(auto const & a : atoms_) {
        if( a == nullptr) { continue; }
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


Residue
str_to_residue(
    String const & s,
    ResidueTypeSet const & rts) {
    Strings spl = split_str_by_delimiter(s, ",");
    ResidueType rtype = rts.get_rtype_by_resname(spl[0]);
    Residue r(rtype, spl[1], std::stoi(spl[2].c_str()), spl[3], spl[4]);
    AtomOPs atoms;
    int i = 5;
    while ( i < spl.size() ) {
        if( spl[i].length() == 1) {
            atoms.push_back( nullptr );
        }
        else {
            AtomOP aop ( new Atom ( str_to_atom(spl[i]) ));
            atoms.push_back(aop);
        }
        i++;
    }
    r.setup_atoms(atoms);
    return r;
}




