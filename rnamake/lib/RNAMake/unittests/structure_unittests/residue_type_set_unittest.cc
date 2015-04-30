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

