//
//  atom.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

//RNAMake Headers
#include "math/xyz_vector.h"
#include "math/numerical.h"
#include "structure/atom.h"

String
Atom::to_str() {
    return name_ + " " + vector_to_str(coords_);
}

String
Atom::to_pdb_str(
    int acount) {

    char buffer [200];
    std::sprintf(buffer, "ATOM %6d  P   C   A   1 %11.3f%8.3f%8.3f  1.00 62.18           P\n", acount, coords_[0], coords_[1], coords_[2]);
    return String(buffer);

}


bool
are_atoms_equal(
        AtomOP const & a1,
        AtomOP const & a2,
        float tol) {

    return are_xyzVector_equal(a1->coords(), a2->coords(), tol) && a1->name() == a2->name();
}


bool
are_atom_vectors_equal(
        AtomOPs const & atoms_1,
        AtomOPs const & atoms_2,
        float tol) {

    if (atoms_1.size() != atoms_2.size()) { return false; }

    int result = 0;
    for (int i = 0; i < atoms_1.size(); i++) {
        result = are_atoms_equal(atoms_1[i], atoms_2[i], tol);
        //std::cout << atoms_1[i]->to_str() << " " << atoms_2[i]->to_str() << std::endl;
        if (!result) { return false; }
    }
    return true;

}