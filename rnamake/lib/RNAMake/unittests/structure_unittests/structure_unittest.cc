//
//  structure_unittest.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 4/30/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "is_equal.hpp"
#include "structure_unittest.h"
#include "unittest.h"
#include "util/file_io.h"
#include "util/settings.h"

namespace unittests {
    
StructureUnittest::StructureUnittest() {
    
    String path = base_dir() + "/rnamake/unittests/resources/motifs/p4p6/p4p6.pdb";
    s_ = std::make_shared<Structure>(path);
    
}

int
StructureUnittest::test_move() {
    auto s2 = Structure(*s_);
    Point p (10.0, 0 ,0);
    s2.move(p);
    AtomOPs new_atoms = s2.atoms();
    AtomOPs org_atoms = s_->atoms();
    
    float dist;
    for(int i = 0; i < org_atoms.size(); i++) {
        if(org_atoms[i].get() == NULL) { continue; }
        dist = 10 - org_atoms[i]->coords().distance(new_atoms[i]->coords());
        if(dist > 0.001) { return 0; }
    }
    
    return 1;
}
    
int
StructureUnittest::test_to_str() {
    ResidueTypeSet rts;
    auto s = s_->to_str();
    auto s_copy = std::make_shared<Structure>(s, rts);
    
    failUnless(are_structures_equal(s_, s_copy, 0), "structures should be equal");
    
    return 1;
}
    

int
StructureUnittest::test_transform() {
    String path = unittest_resource_dir() + "/structure/test_transform.dat";
    Strings lines = get_lines_from_file(path);

    Matrix r = matrix_from_str(lines[0]);
    Point trans = vector_from_str(lines[1]);
    ResidueTypeSet rts;
    auto s2 = Structure(lines[2], rts);

    Transform t(r, trans);
    s_->transform(t);
    
    AtomOPs org_atoms = s2.atoms();
    AtomOPs new_atoms = s_->atoms();
    
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
    

    ResidueOPs res = s_->residues();
    for(auto const & r : res) {
        ResidueOP r_new = s_->get_residue(r->num(), r->chain_id(), r->i_code());
        if(!(r_new->uuid() == r->uuid())) { return 0; }
        r_new = s_->get_residue(r->uuid());
        if(!(r_new->uuid() == r->uuid())) { return 0; }
        
    }
    
    ResidueOP r_new = s_->get_residue(100, "A", "");
    if(r_new.get() != NULL) { return 0; }
    
    return 1;
}

int
StructureUnittest::run() {
    
    if (test_move() == 0)               { std::cout << "test_move failed" << std::endl; }
    if (test_transform() == 0)          { std::cout << "test_transform failed" << std::endl; }
    if (test_get_residue() == 0)        { std::cout << "test_get_residue failed" << std::endl; }
    test_to_str();
    return 0;
}

int
StructureUnittest::run_all() {
    String name = "StructureUnittest";
    typedef int (StructureUnittest::*fptr)();
    std::map<String, fptr> func_map;
    func_map["test_move"              ] = &StructureUnittest::test_move;
    func_map["test_transform"         ] = &StructureUnittest::test_transform;
    func_map["test_get_residue"       ] = &StructureUnittest::test_get_residue;

    int failed = 0;
    for(auto const & kv : func_map) {
        try {
            int result = (this->*kv.second)();
        }
        catch(...) {
            std::cout << name << "::" << kv.first << " returned ERROR!" << std::endl;
            failed += 1;
        }
        
    }
    
    return failed;
}

}







