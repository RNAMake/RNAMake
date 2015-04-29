//
//  atom_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/28/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "atom_unittest.h"

int
AtomUnittest::test_creation() {
    Atom atom("P", Point(0, 1, 2));
    Point p(0, 1, 2);
    if (!are_xyzVector_equal(atom.coords(), p)) {
        return 0;
    }
    return 1;
}

int
AtomUnittest::test_to_pdb_str() {
    Atom atom("H1", Point(1, 2, 3));
    String s = atom.to_pdb_str(1);
    String ref("ATOM      1  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    if (s.compare(ref) != 0) {
        return 1;
    }
    String s2 = atom.to_pdb_str(10);
    String ref2("ATOM     10  P   C   A   1       1.000   2.000   3.000  1.00 62.18           P\n");
    if (s2.compare(ref2) != 0) {
        return 1;
    }
    
    return 1;
}

int
AtomUnittest::test_str_to_atom() {
    Atom atom("H1", Point(1, 2, 3));
    String s = atom.to_str();
    Atom atom2 = str_to_atom(s);
    if(not are_xyzVector_equal(atom.coords(), atom2.coords())) {
        return 0;
    }
    if(not (atom.name().compare(atom2.name()) == 0)) {
        return 0;
    }
    return 1;
}

int
AtomUnittest::test_copy() {
    Atom atom("H1", Point(1, 2, 3));
    Atom atom2 = atom.copy();
    if(not are_xyzVector_equal(atom.coords(), atom2.coords())) {
        return 0;
    }
    if(not (atom.name().compare(atom2.name()) == 0)) {
        return 0;
    }
    return 1;
}

