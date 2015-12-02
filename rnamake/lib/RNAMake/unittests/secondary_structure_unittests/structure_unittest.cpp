//
//  structure_unittest.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 11/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure_unittest.h"
#include "secondary_structure/structure.h"

namespace unittests {
namespace sstruct_unittests  {

void
StructureUnittest::test_creation() {
    auto s = sstruct::Structure("GG+CC", "((+))");
    
    if(s.residues().size() != 4) {
        throw UnittestException("did not get the correct number of residues");
    }
}
 
void
StructureUnittest::test_find_residue() {
    auto ss = sstruct::Structure("AGCU+AGCU", "((((+))))");
    auto r = ss.get_residue(1, "A", "");
    if(r == nullptr) {
        throw UnittestException("did not find a known residue, find_residue not working");
    }
    
    auto r2 = ss.get_residue(r->uuid());
    if(r2 == nullptr) {
        throw UnittestException("did not find a known residue, find_residue not working");
    }
}
    
void
StructureUnittest::test_copy() {
    auto ss = sstruct::Structure("AGCU+AGCU", "((((+))))");
    auto c_ss = sstruct::Structure(ss);
    for(auto const & r : ss.residues()) {
        if(c_ss.get_residue(r->uuid()) == nullptr) {
            throw UnittestException("cannot find residue in copy");
        }
    }
}
    
void
StructureUnittest::test_to_str() {
    auto ss = sstruct::Structure("AGCU+AGCU", "((((+))))");
    auto s = ss.to_str();
    auto c_ss = sstruct::str_to_structure(s);
    for(auto const & r : ss.residues()) {
        if(c_ss.get_residue(r->num(), r->chain_id(), r->i_code()) == nullptr) {
            throw UnittestException("cannot find reisdue in to_str");
        }
    }
}
    
int
StructureUnittest::run() {
    test_creation();
    test_find_residue();
    test_copy();
    test_to_str();
    return 0;
}
    
    
}
}