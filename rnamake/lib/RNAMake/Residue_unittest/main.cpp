//
//  main.cpp
//  Residue_unittest
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "residue.h"
#include "residue_type_set.h"

int
test_bead_creation() {
    Bead b ( Point(1, 2, 3), BASE);
    return 1;
}

int
test_str_to_residue() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    while ( input.good() ) {
        getline(input, line);
        if( line.length() < 10 ) { break; }
        Residue r = str_to_residue(line, rts);
    }
    
    return 1;
}

int
test_get_atom() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    getline(input, line);
    ResidueTypeSet rts;
    Residue r = str_to_residue(line, rts);
    String name = "C1'";
    AtomOP a = r.get_atom(name);
    std::cout << a->name() << std::endl;
    return 1;
}


int main(int argc, const char * argv[]) {
    if (test_bead_creation() == 0)   { std::cout << "test_bead_creation failed" << std::endl; }
    //if (test_str_to_residue() == 0)  { std::cout << "test_str_to_residue failed" << std::endl; }
    if (test_get_atom() == 0)        { std::cout << "test_get_atom failed" << std::endl; }
    
    return 0;
}
