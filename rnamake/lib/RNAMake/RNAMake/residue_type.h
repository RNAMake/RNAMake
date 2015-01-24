//
//  residue_type.h
//  RNAMake
//
//  Created by Joseph Yesselman on 1/24/15.
//  Copyright (c) 2015 Joseph Yesselman. All rights reserved.
//

#ifndef __RNAMake__residue_type__
#define __RNAMake__residue_type__

#include <stdio.h>
#include "types.h"
#include "atom.h"

class ResidueType {
public:
    ResidueType() {}
    ResidueType(
        String const &,
        StringIntMap const &);
    
    ~ResidueType() {}

public:
    StringOP
    get_correct_atom_name(
        Atom const &);

private:
    String name_;
    StringIntMap atom_map_;
    Strings alt_names_;

};


#endif /* defined(__RNAMake__residue_type__) */
