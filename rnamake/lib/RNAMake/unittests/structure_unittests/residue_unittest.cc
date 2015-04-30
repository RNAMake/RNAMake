//
//  residue_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/file_io.h"

#include "unittest.h"
#include "residue_unittest.h"



int
ResidueUnittest::test_bead_creation() {
    Bead b ( Point(1, 2, 3), BASE);
    return 1;
}

int
ResidueUnittest::test_str_to_residue() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);
    ResidueTypeSet rts;
    for(auto const & line : lines) {
        Residue r = str_to_residue(line, rts);
    }
    
    return 1;
}

int
ResidueUnittest::test_get_atom() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);
    ResidueTypeSet rts;
    Residue r = str_to_residue(lines[0], rts);
    String name = "C1'";
    AtomOP a = r.get_atom(name);
    return 1;
}

int
ResidueUnittest::test_connected_to() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    Residue r1 = str_to_residue(lines[0], rts);
    Residue r2 = str_to_residue(lines[1], rts);
    Residue r3 = str_to_residue(lines[2], rts);
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
ResidueUnittest::test_get_beads() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    Residue r = str_to_residue(lines[0], rts);
    Beads beads = r.get_beads();
    return 1;
}

int
ResidueUnittest::test_copy() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    Residue r = str_to_residue(lines[0], rts);
    Residue r2 = r.copy();
    
    return 1;
}

int
ResidueUnittest::test_to_str() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    Residue r = str_to_residue(lines[0], rts);
    String s = r.to_str();
    Residue r2 = str_to_residue(s, rts);
    String name = "N1";
    AtomOP a = r2.get_atom(name);
    return 1;
}

int
ResidueUnittest::test_equals() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    Residue r1 = str_to_residue(lines[0], rts);
    Residue r2 = str_to_residue(lines[1], rts);

    if(!(r1 == r1)) { return 0; }
    if(r1 == r2) { return 0; }
    return 1;
}


int
ResidueUnittest::test_memory_management() {
    String file = "test_str_to_residue.dat";
    String line;
    std::ifstream input;
    input.open(file);
    ResidueTypeSet rts;
    getline(input, line);
    input.close();
    for(int i = 0; i < 100000; i++) {
        Residue r = str_to_residue(line, rts);
    }
    return 1;
}


int
ResidueUnittest::run() {
    if (test_bead_creation() == 0)   { std::cout << "test_bead_creation failed" << std::endl; }
    if (test_str_to_residue() == 0)  { std::cout << "test_str_to_residue failed" << std::endl; }
    if (test_get_atom() == 0)        { std::cout << "test_get_atom failed" << std::endl; }
    if (test_connected_to() == 0)    { std::cout << "test_connected_to failed" << std::endl; }
    if (test_get_beads() == 0)       { std::cout << "test_get_beads failed" << std::endl; }
    if (test_copy() == 0)            { std::cout << "test_copy failed" << std::endl; }
    if (test_to_str() == 0)          { std::cout << "test_to_str failed" << std::endl; }
    if (test_equals() == 0)          { std::cout << "test_equals failed" << std::endl; }
    //test_memory_management();
    
    
    
    return 0;
}