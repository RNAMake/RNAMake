//
//  main.cpp
//  Residue_unittest
//
//  Created by Joseph Yesselman on 1/25/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "uuid.h"
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
    return 1;
}

int
test_connected_to() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Residue r1 = str_to_residue(line, rts);
    getline(input, line);
    Residue r2 = str_to_residue(line, rts);
    getline(input, line);
    Residue r3 = str_to_residue(line, rts);
    if (r1.connected_to(r2, 3.0)  != 1) {
        return 0;
    }
    if (r1.connected_to(r3, 3.0)  != 0) {
        return 0;
    }
    if (r3.connected_to(r2, 3.0)  != -1) {
        return 0;
    }
    
    return 1;
}

int
test_get_beads() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    getline(input, line);
    ResidueTypeSet rts;
    Residue r = str_to_residue(line, rts);
    Beads beads = r.get_beads();
    return 1;
}

int
test_copy() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    getline(input, line);
    ResidueTypeSet rts;
    Residue r = str_to_residue(line, rts);
    Residue r2 = r.copy();
    
    return 1;
}

int
test_to_str() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    getline(input, line);
    ResidueTypeSet rts;
    Residue r = str_to_residue(line, rts);
    String s = r.to_str();
    Residue r2 = str_to_residue(s, rts);
    String name = "N1";
    AtomOP a = r2.get_atom(name);
    //int acount = 1;
    //std::cout << acount << std::endl;
    return 1;
}

int
test_uuid() {
    Uuid uuid1, uuid2;
    if (uuid1 == uuid2) { return 0; }
    if (!(uuid1 == uuid1)) { return 0; }
    return 1;
}

int
test_equals() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    Residue r1 = str_to_residue(line, rts);
    getline(input, line);
    Residue r2 = str_to_residue(line, rts);
    std::cout << (r1 == r1) << std::endl;
    std::cout << (r1 == r2) << std::endl;
    return 1;
}


int main(int argc, const char * argv[]) {
    if (test_uuid() == 0)            { std::cout << "test_uuid failed" << std::endl; }
    if (test_bead_creation() == 0)   { std::cout << "test_bead_creation failed" << std::endl; }
    if (test_str_to_residue() == 0)  { std::cout << "test_str_to_residue failed" << std::endl; }
    if (test_get_atom() == 0)        { std::cout << "test_get_atom failed" << std::endl; }
    if (test_connected_to() == 0)    { std::cout << "test_connected_to failed" << std::endl; }
    if (test_get_beads() == 0)       { std::cout << "test_get_beads failed" << std::endl; }
    if (test_copy() == 0)            { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)          { std::cout << "test_to_str failed" << std::endl; }
    if (test_equals() == 0)          { std::cout << "test_equals failed" << std::endl; }

    return 0;
}
