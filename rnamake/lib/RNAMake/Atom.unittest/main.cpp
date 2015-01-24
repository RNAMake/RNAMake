//
//  main.cpp
//  Atom.unittest
//
//  Created by Joseph Yesselman on 12/23/14.
//  Copyright (c) 2014 Joseph Yesselman. All rights reserved.
//

#include "gtest/gtest.h"
#include <iostream>
#include "xyzVector.h"
#include "atom.h"
#include "numerical.h"


int
test_creation() {
    Atom atom("P", Point(0, 1, 2));
    Point p(0, 1, 2);
    if (!are_xyzVector_equal(atom.coords(), p)) {
        return 0;
    }
    return 1;
}

int
test_to_pdb_str() {
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
test_str_to_atom() {
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
test_copy() {
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


int main(int argc, const char * argv[]) {
    if (test_creation() == 0)    {  std::cout << "test_creation failed" << std::endl; }
    if (test_to_pdb_str() == 0)  {  std::cout << "test_to_pdb_str failed" << std::endl; }
    if (test_str_to_atom() == 0) {  std::cout << "test_str_to_atom failed" << std::endl; }
    if (test_copy() == 0)        {  std::cout << "test_copy failed" << std::endl; }

}








