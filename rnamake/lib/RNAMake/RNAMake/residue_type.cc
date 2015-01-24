//
//  residue_type.cc
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#include "residue_type.h"
#include "types.h"
#include "atom.h"

ResidueType::ResidueType(
    String const & name,
    StringIntMap const & atom_map):
    name_ (name),
    atom_map_ ( atom_map) {
    
    alt_names_ = Strings();
    alt_names_.push_back(name_.substr(1));
        
}

StringOP
ResidueType::get_correct_atom_name(
    Atom const & a) {

    if( atom_map.find(a.name()) != atom_map_.end()) {
        return StringOP(atom_map_[a.name()]);
    }
    else {
        return NULL;
    }
}











