//
//  is_equal.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "is_equal.hpp"
#include "math/numerical.h"


bool
are_atoms_equal(
    AtomOP const & a1,
    AtomOP const & a2,
    float tol) {
    
    return are_xyzVector_equal(a1->coords(), a2->coords(), tol) && a1->name() == a2->name();
}


bool
are_atom_vectors_equal(
    AtomOPs const & atoms_1,
    AtomOPs const & atoms_2,
    float tol) {
    
    if(atoms_1.size() != atoms_2.size()) { return false; }
    
    int result = 0;
    for(int i = 0; i < atoms_1.size(); i++) {
        result = are_atoms_equal(atoms_1[i], atoms_2[i], tol);
        //std::cout << atoms_1[i]->to_str() << " " << atoms_2[i]->to_str() << std::endl;
        if(!result) { return false; }
    }
    return true;
    
}

bool
are_residues_equal(
    ResidueOP const & r1,
    ResidueOP const & r2,
    int check_uuids) {
    
    if(r1->name() != r2->name()) { return false; }
    if(check_uuids && r1->uuid() != r2->uuid()) { return false; }
    
    int i = -1;
    bool result;
    for(auto const & a : *r1) {
        i++;
        if(a == nullptr && r2->has_atom(i))   { return 0; }
        if(!r2->has_atom(i) && a != nullptr)  { return 0;}
        if(a == nullptr and !r2->has_atom(i)) { continue; }
        
        result = are_atoms_equal(a, r2->get_atom(i));
        if(!result) { return false; }
    }
    
    return true;
    
}

bool
are_chains_equal(
    ChainOP const & c1,
    ChainOP const & c2,
    int check_uuids) {
    
    if(c1->length() != c2->length()) { return false; }

    auto result = 0;
    for(int i = 0; i < c1->length(); i++) {
        result = are_residues_equal(c1->get_residue(i), c2->get_residue(i), check_uuids);
        if(!result) { return false; }
    }
    return true;
}


bool
are_structures_equal(
    StructureOP const & s1,
    StructureOP const & s2,
    int check_uuids) {

    if(s1->num_chains() != s2->num_chains()) { return false; }
    for(int i = 0; i < s1->num_chains(); i++) {
        auto result = are_chains_equal(s1->get_chain(i), s2->get_chain(i), check_uuids);
        if(!result) { return false; }
    }

    return true;
    
}
