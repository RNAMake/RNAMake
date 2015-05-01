//
//  structure_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "structure_unittest.h"
#include "unittest.h"
#include "util/file_io.h"

StructureUnittest::StructureUnittest() {
    
    String path = unittest_resource_dir() + "/structure/test_str_to_structure.dat";
    Strings lines = get_lines_from_file(path);
    
    ResidueTypeSet rts;
    s_ = str_to_structure(lines[0], rts);
}

int
StructureUnittest::test_build_chains() {
    Structure s2 = s_.copy();
    ResidueOPs res = s_.residues();
    s_._build_chains(res);
    for(int i = 0; i < s_.residues().size(); i++ ){
        if(s_.residues()[i]->num() != s2.residues()[i]->num()) { return 0; }
    }
    
    return 1;

}

int
StructureUnittest::test_move() {
    Structure s2 = s_.copy();
    Point p (10.0, 0 ,0);
    s2.move(p);
    AtomOPs new_atoms = s2.atoms();
    AtomOPs org_atoms = s_.atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = 10 - org_atoms[i]->coords().distance(new_atoms[i]->coords());
        if(dist > 0.001) { return 0; }
    }
    
    return 1;
    
}

int
StructureUnittest::test_transform() {
    String path = unittest_resource_dir() + "/structure/test_transform.dat";
    Strings lines = get_lines_from_file(path);

    Matrix r = matrix_from_str(lines[0]);
    Point trans = vector_from_str(lines[1]);
    ResidueTypeSet rts;
    Structure s2 = str_to_structure(lines[2], rts);

    Transform t(r, trans);
    s_.transform(t);
    
    AtomOPs org_atoms = s2.atoms();
    AtomOPs new_atoms = s_.atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = org_atoms[i]->coords().distance(new_atoms[i]->coords());
        if(dist > 0.001) { return 0; }
    }
    
    return 1;
}

int
StructureUnittest::test_get_residue() {
    

    ResidueOPs res = s_.residues();
    for(auto const & r : res) {
        ResidueOP r_new = s_.get_residue(r->num(), r->chain_id(), r->i_code());
        if(!(r_new->uuid() == r->uuid())) { return 0; }
        r_new = s_.get_residue(r->uuid());
        if(!(r_new->uuid() == r->uuid())) { return 0; }
        
    }
    
    ResidueOP r_new = s_.get_residue(100, "A", "");
    if(r_new.get() != NULL) { return 0; }
    
    return 1;
}


int
StructureUnittest::run() {
    if (test_build_chains() == 0)       { std::cout << "test_build_chains failed" << std::endl; }
    if (test_move() == 0)               { std::cout << "test_move failed" << std::endl; }
    if (test_transform() == 0)          { std::cout << "test_transform failed" << std::endl; }
    if (test_get_residue() == 0)        { std::cout << "test_get_residue failed" << std::endl; }

    return 0;
}







