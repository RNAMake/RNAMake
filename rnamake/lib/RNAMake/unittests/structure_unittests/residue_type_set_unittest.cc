//
//  residue_type_set_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/29/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue_type_set_unittest.h"


int
ResidueTypeSetUnittest::test_creation_residue_type() {
    String name("GUA");
    StringIntMap atom_map;
    atom_map["P"] = 0;
    
    ResidueType rt (name, atom_map);
    return 1;
}

int
ResidueTypeSetUnittest::test_match_name() {
    String name("GUA");
    StringIntMap atom_map;
    atom_map["P"] = 0;
    
    ResidueType rt (name, atom_map);
    String name1 = "GUA", name2 = "G", name3 = "rG", name4 = "rC";
    int result = rt.match_name(name1);
    if (!result) { return 0; }
    result = rt.match_name(name2);
    if (!result) { return 0; }
    result = rt.match_name(name3);
    if (!result) { return 0; }
    result = rt.match_name(name4);
    if (result) { return 0; }
    return 1;
}

int
ResidueTypeSetUnittest::test_creation() {
    ResidueTypeSet rts;
    return 1;
}

int
ResidueTypeSetUnittest::test_get_rtype_by_resname() {
    ResidueTypeSet rts;
    String name1 = "GUA", name2 = "@%$";
    ResidueType rt = rts.get_rtype_by_resname(name1);
    if(rt.short_name().compare("G") != 0 ) {
        return 0;
    }
    try {
        ResidueType rt2= rts.get_rtype_by_resname(name2);
        return 0;
    }
    catch(String e) { }
    
    return 1;
}

int
ResidueTypeSetUnittest::test_atom_pos_by_name() {
    ResidueTypeSet rts;
    ResidueType rtype = rts.get_rtype_by_resname("GUA");
    String name = "P";
    rtype.atom_pos_by_name(name);
    return 1;
}

int
ResidueTypeSetUnittest::run() {
    if (test_creation_residue_type() == 0) { std::cout << "test_creation_residue_type failed" << std::endl;  }
    if (test_creation() == 0 )             { std::cout << "test_creation failed" << std::endl; }
    if (test_match_name() == 0 )           { std::cout << "test_match_name failed" << std::endl; }
    if (test_get_rtype_by_resname() == 0 ) { std::cout << "test_get_rtype_by_resname failed" << std::endl; }
    if (test_atom_pos_by_name() == 0 )     { std::cout << "test_atom_pos_by_name failed" << std::endl; }
    
    return 1;
}

int
ResidueTypeSetUnittest::run_all() {
    String name = "ResidueTypeSetUnittest";
    typedef int (ResidueTypeSetUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_creation_residue_type" ] = &ResidueTypeSetUnittest::test_creation_residue_type;
    func_map["test_match_name"            ] = &ResidueTypeSetUnittest::test_match_name;
    func_map["test_creation"              ] = &ResidueTypeSetUnittest::test_creation;
    func_map["test_get_rtype_by_resname"  ] = &ResidueTypeSetUnittest::test_get_rtype_by_resname;
    func_map["test_atom_pos_by_name"      ] = &ResidueTypeSetUnittest::test_atom_pos_by_name;

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


















