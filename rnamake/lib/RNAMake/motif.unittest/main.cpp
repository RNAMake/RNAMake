//
//  main.cpp
//  motif.unittest
//
//  Created by Joseph Yesselman on 1/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "motif.h"
#include "residue_type_set.h"

Motif
get_test_motif() {
    String file = "test_str_to_motif.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Motif m(line, rts);
    return m;
}

int
test_str_to_motif() {
    String file = "test_str_to_motif.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Motif m(line, rts);
    
    return 1;
}

int
test_copy() {
    Motif m = get_test_motif();
    Motif mcopy = m.copy();
    Point p(10,20,20);
    m.move(p);
    m.to_pdb("original.pdb");
    mcopy.to_pdb("copy.pdb");
    return 1;
}

int
test_to_str() {
    Motif m = get_test_motif();
    ResidueTypeSet rts;
    Motif m2 = Motif(m.to_str(), rts);
    return 1;
}

int main(int argc, const char * argv[]) {
    if (test_str_to_motif() == 0)     { std::cout << "test_str_to_motif failed" << std::endl; }
    if (test_copy() == 0)             { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)           { std::cout << "test_to_str failed" << std::endl; }
    return 0;
}
