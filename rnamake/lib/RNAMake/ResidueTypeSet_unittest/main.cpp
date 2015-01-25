//
//  main.cpp
//  ResidueTypeSet_unittest
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include <iostream>
#include "residue_type_set.h"
#include "residue_type.h"
#include "types.h"

int
test_creation_residue_type() {
    String name("GUA");
    StringIntMap atom_map;
    atom_map["P"] = 0;
    
    ResidueType rt (name, atom_map);
    return 1;
}

int
test_match_name() {
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
test_creation() {
    ResidueTypeSet rts;
    return 1;
}

int
test_get_rtype_by_resname() {
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


int main(int argc, const char * argv[]) {
    if (test_creation_residue_type() == 0) { std::cout << "test_creation_residue_type failed" << std::endl;  }
    if (test_creation() == 0 )             { std::cout << "test_creation failed" << std::endl; }
    if (test_match_name() == 0 )           { std::cout << "test_match_name failed" << std::endl; }
    if (test_get_rtype_by_resname() == 0 ) { std::cout << "test_get_rtype_by_resname failed" << std::endl; }
    
    return 0;
}
