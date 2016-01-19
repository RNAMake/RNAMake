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

namespace unittests {

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
        if(line.size() < 10) { break; }
        auto r = Residue(line, rts);
    }
    
    return 1;
}

int
ResidueUnittest::test_get_atom() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);
    ResidueTypeSet rts;
    auto r = Residue(lines[0], rts);
    String name = "C1'";
    AtomOP a = r.get_atom(name);
    return 1;
}

int
ResidueUnittest::test_connected_to() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    auto r1 = Residue(lines[0], rts);
    auto r2 = Residue(lines[1], rts);
    auto r3 = Residue(lines[2], rts);
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
    auto r = Residue(lines[0], rts);
    Beads beads = r.get_beads();
    return 1;
}

int
ResidueUnittest::test_copy() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    auto r  = Residue(lines[0], rts);
    auto r2 = Residue(r);
    
    return 1;
}

int
ResidueUnittest::test_to_str() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    auto r = Residue(lines[0], rts);
    String s = r.to_str();
    auto r2 = Residue(s, rts);
    String name = "N1";
    AtomOP a = r2.get_atom(name);
    return 1;
}

int
ResidueUnittest::test_equals() {
    String path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    Strings lines = get_lines_from_file(path);

    ResidueTypeSet rts;
    auto r1 = Residue(lines[0], rts);
    auto r2 = Residue(lines[1], rts);

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
        auto r = Residue(line, rts);
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

int
ResidueUnittest::run_all() {
    String name = "ResidueUnittest";
    typedef int (ResidueUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_bead_creation"  ] = &ResidueUnittest::test_bead_creation;
    func_map["test_str_to_residue" ] = &ResidueUnittest::test_str_to_residue;
    func_map["test_get_atom"       ] = &ResidueUnittest::test_get_atom;
    func_map["test_connected_to"   ] = &ResidueUnittest::test_connected_to;
    func_map["test_get_beads"      ] = &ResidueUnittest::test_get_beads;
    func_map["test_copy"           ] = &ResidueUnittest::test_copy;
    func_map["test_to_str"         ] = &ResidueUnittest::test_to_str;
    func_map["test_equals"         ] = &ResidueUnittest::test_equals;
    
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
            if(result == 0) {
                std::cout << name << "::" << kv.first << " FAILED!" << std::endl;
            }
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }
    
    return 0;
}


}













