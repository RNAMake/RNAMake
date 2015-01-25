//
//  residue.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue.h"

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
