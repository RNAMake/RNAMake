//
//  residue_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "util/file_io.h"

#include "is_equal.hpp"
#include "unittest.h"
#include "residue_unittest.h"

namespace unittests {

ResidueUnittest::ResidueUnittest():
rts_(ResidueTypeSet()),
residues_(ResidueOPs()){
    auto path = unittest_resource_dir() + "residue/test_str_to_residue.dat";
    auto lines = get_lines_from_file(path);
    for(auto const & line : lines) {
        if(line.size() < 10) { break; }
        auto r = std::make_shared<Residue>(line, rts_);
        residues_.push_back(r);
    }
}
    
void
ResidueUnittest::test_bead_creation() {
    auto b =  Bead( Point(1, 2, 3), BASE);
}


void
ResidueUnittest::test_get_atom() {
    auto r = residues_[0];
    auto name = String("C1'");
    auto a = r->get_atom(name);
}

void
ResidueUnittest::test_connected_to() {
    auto r1 = residues_[0];
    auto r2 = residues_[1];
    auto r3 = residues_[2];
    
    failUnless(r1->connected_to(*r2) == 1,  "r1 should be connected to r2");
    failUnless(r1->connected_to(*r3) == 0,  "r1 should not be connected to r3");
    failUnless(r3->connected_to(*r2) == -1, "r3 should be connected to r2");

}

void
ResidueUnittest::test_get_beads() {
    auto r = residues_[0];
    Beads beads = r->get_beads();
}

void
ResidueUnittest::test_copy() {
    auto r  = residues_[0];
    auto r2 = std::make_shared<Residue>(*r);
    failUnless(are_residues_equal(r, r2, 0), "residues should be equal but arent");
}

void
ResidueUnittest::test_to_str() {
    auto r = residues_[0];
    auto s = r->to_str();
    auto r2 = std::make_shared<Residue>(s, rts_);

    failUnless(are_residues_equal(r, r2, 0), "residues should be equal but arent");
}

void
ResidueUnittest::test_equals() {
    auto r1 = residues_[0];
    auto r2 = residues_[1];
    auto r1_copy = std::make_shared<Residue>(*r1);
    
    failUnless(!(*r1 == *r2), "r1 and r2 are not equal");
    failUnless(*r1 == *r1_copy, "r1 and r1_copy are equal");
}

void
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
}

    
int
ResidueUnittest::run() {
    test_bead_creation();
    test_get_atom();
    test_connected_to();
    test_get_beads();
    test_copy();
    test_to_str();
    test_equals();
    //test_memory_management();
    
    return 0;
}

int
ResidueUnittest::run_all() {
    String name = "ResidueUnittest";
    typedef void (ResidueUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_bead_creation"  ] = &ResidueUnittest::test_bead_creation;
    func_map["test_get_atom"       ] = &ResidueUnittest::test_get_atom;
    func_map["test_connected_to"   ] = &ResidueUnittest::test_connected_to;
    func_map["test_get_beads"      ] = &ResidueUnittest::test_get_beads;
    func_map["test_copy"           ] = &ResidueUnittest::test_copy;
    func_map["test_to_str"         ] = &ResidueUnittest::test_to_str;
    func_map["test_equals"         ] = &ResidueUnittest::test_equals;
    
    for(auto const & kv : func_map) {
        try {
            (this->*kv.second)();
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
        }
        
    }
    
    return 0;
}


}













