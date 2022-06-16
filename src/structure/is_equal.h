//
//  is_equal.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 4/17/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#ifndef is_equal_hpp
#define is_equal_hpp

#include <stdio.h>

//#include "structure/structure.h"
#include "structure/atom.h"
#include "structure/residue.h"
#include <math/vector_3.hpp>

namespace structure {

bool are_atoms_equal(AtomOP const &, AtomOP const &, float tol = 0.001);

bool are_atom_vectors_equal(AtomOPs const &, AtomOPs const &,
                            float tol = 0.001);

bool are_residues_equal(ResidueOP const &r1, ResidueOP const &r2,
                        int check_uuids = 1);

// bool
// are_chains_equal(
//         ChainOP const & c1,
//         ChainOP const & c2,
//         int check_uuids = 1);

// bool
// are_structures_equal(
//         StructureOP const & s1,
//         StructureOP const & s2,
//         int check_uuids = 1);

} // namespace structure

#endif /* is_equal_hpp */
