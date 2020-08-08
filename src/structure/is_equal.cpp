//
//  is_equal.cpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include "structure/is_equal.h"
#include "math/numerical.h"

namespace structure {

bool
are_atoms_equal(
        AtomOP const & a1,
        AtomOP const & a2,
        float tol) {

    return math::are_xyzVector_equal(a1->coords(), a2->coords(), tol) &&
           a1->name() == a2->name();
}


bool
are_atom_vectors_equal(
        AtomOPs const & atoms_1,
        AtomOPs const & atoms_2,
        float tol) {

    if (atoms_1.size() != atoms_2.size()) { return false; }

    for (int i = 0; i < atoms_1.size(); i++) {
        int result = are_atoms_equal(atoms_1[i], atoms_2[i], tol);
        //std::cout << atoms_1[i]->to_str() << " " << atoms_2[i]->to_str() << std::endl;
        if (!result) { return false; }
    }
    return true;

}

bool
are_residues_equal(
        ResidueOP const & r1,
        ResidueOP const & r2,
        int check_uuids) {

    if (r1->name() != r2->name()) { return false; }
    if (check_uuids && r1->uuid() != r2->uuid()) { return false; }

    int i = -1;
    auto r2_atoms = r2->atoms();
    bool result;
    for (auto const & a : r1->atoms()) {
        i++;
        if (a == nullptr && r2_atoms[i] != nullptr) { return false; }
        else if (a != nullptr && r2_atoms[i] == nullptr) { return false; }
        else if (a == nullptr && r2_atoms[i] == nullptr) { continue; }

        result = are_atoms_equal(a, r2_atoms[i]);
        if (!result) { return false; }

    }

    return true;

}


bool
are_chains_equal(
        ChainOP const & c1,
        ChainOP const & c2,
        int check_uuids) {

    if (c1->length() != c2->length()) { return false; }

    auto c1_res = c1->residues();
    auto c2_res = c2->residues();
    auto result = 0;
    for (int i = 0; i < c1->length(); i++) {
        result = are_residues_equal(c1_res[i], c2_res[i], check_uuids);
        if (!result) { return false; }
    }
    return true;
}

bool
are_structures_equal(
        StructureOP const & s1,
        StructureOP const & s2,
        int check_uuids) {
    
    const auto& chains1 = s1->chains();
    const auto& chains2 = s2->chains();
    
    if(chains1.size() != chains2.size()) {
        return false;
    }

    auto it_s1 = chains1.cbegin();
    auto it_s2 = chains2.cbegin();
    const auto stop = chains1.cend();
    for( ; it_s1 != stop; ++it_s1, ++it_s2) {
        if(!are_chains_equal(
                    *it_s1,
                    *it_s2
                    )) {
            return false;
        }
    }

    return true;

}

}
